using AdaptiveOpticsProperHIL
using AdaptiveOpticsSim
import AMDGPU
import CUDA
using SPIDERSSim

using AdaptiveOpticsSim.AlgorithmGraphs

const Backends = AdaptiveOpticsSim.Backends
const AOS_ROOT = dirname(dirname(pathof(AdaptiveOpticsSim)))

include(joinpath(
    AOS_ROOT,
    "examples",
    "integrations",
    "pyrtc",
    "pyrtc_shared_memory.jl",
))
using .PyRTCSharedMemory

const STREAM_NAMES = (
    "spidersPupilOpd",
    "spidersLlowfs",
    "spidersSccFrame",
    "spidersSccDifference",
)

function load_backend(name::Symbol)
    if name === :cpu
        return Backends.HostComputeDevice()
    elseif name === :amdgpu
        AMDGPU.functional() || error("AMDGPU is not functional on this host")
        AMDGPU.allowscalar(false)
        probe = AMDGPU.ROCArray(zeros(Float32, 1))
        return Backends.compute_device(probe)
    elseif name === :cuda
        CUDA.functional() || error("CUDA is not functional on this host")
        CUDA.allowscalar(false)
        probe = CUDA.CuArray(zeros(Float32, 1))
        return Backends.compute_device(probe)
    end
    throw(ArgumentError("backend must be cpu, amdgpu, or cuda"))
end

function device_copy(target, source::Matrix{T}) where {T}
    destination = Backends.allocate_device_array(target, T, size(source)...)
    copyto!(destination, source)
    return destination
end

function graph_definition(configuration, pupil_opd, pupil_amplitude, name)
    return algorithm_graph(
        (proper_propagation_node(:proper, configuration),);
        name,
        inputs=(
            graph_input(:pupil_opd, :proper => :pupil_opd, pupil_opd),
            graph_input(
                :pupil_amplitude,
                :proper => :pupil_amplitude,
                pupil_amplitude,
            ),
        ),
        outputs=(graph_output(:intensity, :proper => :output),),
    )
end

function moving_opd!(opd::Matrix{T}, frame::Int) where {T}
    phase = T(frame) * T(0.035)
    rows, columns = size(opd)
    row_center = T(rows + 1) / T(2)
    column_center = T(columns + 1) / T(2)
    two_pi = T(2pi)
    @inbounds for column in axes(opd, 2)
        x = (T(column) - column_center) / T(columns)
        for row in axes(opd, 1)
            y = (T(row) - row_center) / T(rows)
            opd[row, column] = T(2.5e-7) * (
                T(0.55) * sin(two_pi * (x + phase)) +
                T(0.30) * sin(two_pi * (T(2) * y - T(0.7) * phase)) +
                T(0.15) * sin(two_pi * (x + y + T(1.3) * phase))
            )
        end
    end
    return opd
end

function viewer_executable()
    if haskey(ENV, "PYRTC_VIEWER")
        executable = abspath(ENV["PYRTC_VIEWER"])
        isfile(executable) || error("PYRTC_VIEWER does not name a file")
        return executable
    end
    python = get(ENV, "PYRTC_PYTHON", "")
    isempty(python) && error(
        "set PYRTC_VIEWER to pyrtc-view or PYRTC_PYTHON to the pyRTC " *
        "environment's Python interpreter",
    )
    executable = joinpath(dirname(abspath(python)), "pyrtc-view")
    isfile(executable) || error("the selected environment has no pyrtc-view")
    return executable
end

function close_stream_noexcept!(stream)
    isnothing(stream) && return nothing
    try
        unlink!(stream)
    catch
    end
    try
        close(stream)
    catch
    end
    return nothing
end

function run_demo(
    backend_name::Symbol;
    duration::Real=120,
    frame_rate::Real=20,
)
    duration > 0 || throw(ArgumentError("duration must be positive"))
    frame_rate > 0 || throw(ArgumentError("frame_rate must be positive"))

    target = load_backend(backend_name)
    profile = provisional_spiders_h_regular_profile(Float32)
    llowfs_configuration = provisional_spiders_llowfs_configuration(profile)
    scc_configuration = provisional_spiders_chopped_scc_configuration(profile)
    pupil_shape = (
        llowfs_configuration.pupil_resolution,
        llowfs_configuration.pupil_resolution,
    )
    pupil_opd_host = zeros(Float32, pupil_shape)
    pupil_amplitude_host = ones(Float32, pupil_shape)
    pupil_opd = device_copy(target, pupil_opd_host)
    pupil_amplitude = device_copy(target, pupil_amplitude_host)
    llowfs_graph = prepare_algorithm_graph(
        graph_definition(
            llowfs_configuration,
            pupil_opd,
            pupil_amplitude,
            :spiders_llowfs_viewer,
        );
        target,
        execution=StreamGraphExecution(),
    )
    scc_graph = prepare_algorithm_graph(
        graph_definition(
            scc_configuration,
            pupil_opd,
            pupil_amplitude,
            :spiders_scc_viewer,
        );
        target,
        execution=StreamGraphExecution(),
    )

    llowfs = zeros(Float32, profile.llowfs_output_shape)
    scc_frame = zeros(Float32, profile.scc_output_shape)
    scc_fringed = similar(scc_frame)
    scc_unfringed = similar(scc_frame)
    scc_difference = zeros(Float32, profile.scc_output_shape)
    have_fringed = false
    have_unfringed = false
    pupil_stream = nothing
    llowfs_stream = nothing
    scc_stream = nothing
    difference_stream = nothing
    viewer = nothing
    try
        pupil_stream = create_stream(STREAM_NAMES[1], Float32, pupil_shape)
        llowfs_stream = create_stream(
            STREAM_NAMES[2],
            Float32,
            profile.llowfs_output_shape,
        )
        scc_stream = create_stream(
            STREAM_NAMES[3],
            Float32,
            profile.scc_output_shape,
        )
        difference_stream = create_stream(
            STREAM_NAMES[4],
            Float32,
            profile.scc_output_shape,
        )
        command = Cmd(String[
            viewer_executable(),
            STREAM_NAMES...,
            "--geometry",
            "2x2",
            "--fps",
            string(max(1, floor(Int, frame_rate / 2))),
            "--pixel-scale",
            "2",
        ])
        if !haskey(ENV, "QT_QPA_PLATFORM") &&
                !haskey(ENV, "DISPLAY") && haskey(ENV, "WAYLAND_DISPLAY")
            command = addenv(command, "QT_QPA_PLATFORM" => "wayland")
        end
        viewer = run(command; wait=false)

        started = time()
        next_frame = started
        frame = 0
        fringed_count = 0
        unfringed_count = 0
        while time() - started < duration && process_running(viewer)
            frame += 1
            moving_opd!(pupil_opd_host, frame)
            copyto!(pupil_opd, pupil_opd_host)
            step_graph!(llowfs_graph)
            step_graph!(scc_graph)
            copyto!(llowfs, graph_output(llowfs_graph, :intensity))
            copyto!(scc_frame, graph_output(scc_graph, :intensity))

            status = spiders_chopper_frame_status(scc_graph)
            phase = spiders_chopper_phase(status)
            if phase === SpidersFringed
                fringed_count += 1
                copyto!(scc_fringed, scc_frame)
                have_fringed = true
            else
                unfringed_count += 1
                copyto!(scc_unfringed, scc_frame)
                have_unfringed = true
            end
            if have_fringed && have_unfringed
                @. scc_difference = scc_fringed - scc_unfringed
            end

            publish!(pupil_stream, pupil_opd_host)
            publish!(llowfs_stream, llowfs)
            publish!(scc_stream, scc_frame)
            publish!(difference_stream, scc_difference)

            if frame == 1 || frame % max(1, round(Int, frame_rate)) == 0
                println(
                    "frame ", frame,
                    ", SCC phase ", phase,
                    ", LLOWFS range ", extrema(llowfs),
                    ", SCC difference range ", extrema(scc_difference),
                )
            end
            next_frame += inv(Float64(frame_rate))
            delay = next_frame - time()
            delay > 0 && sleep(delay)
        end
        println(
            "SPIDERS viewer completed ", frame, " frames; SCC counts ",
            "fringed=", fringed_count, ", unfringed=", unfringed_count,
        )
    finally
        if !isnothing(viewer) && process_running(viewer)
            kill(viewer)
            wait(viewer)
        end
        close_stream_noexcept!(difference_stream)
        close_stream_noexcept!(scc_stream)
        close_stream_noexcept!(llowfs_stream)
        close_stream_noexcept!(pupil_stream)
    end
    return nothing
end

function main(args)
    backend = Symbol(isempty(args) ? "cpu" : args[1])
    duration = length(args) >= 2 ? parse(Float64, args[2]) : 120.0
    frame_rate = length(args) >= 3 ? parse(Float64, args[3]) : 20.0
    run_demo(backend; duration, frame_rate)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
