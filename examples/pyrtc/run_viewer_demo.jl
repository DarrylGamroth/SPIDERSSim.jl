using AdaptiveOpticsProperHIL
using AdaptiveOpticsSim
import AMDGPU
import CUDA
using SPIDERSSim

using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Optics: Source, photon_irradiance

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
    "spidersPupilAtmosphereOpd",
    "spidersLlowfsPhotonRate",
    "spidersSccUnfringedPhotonRate",
    "spidersSccDifferencePhotonRate",
    "spidersDmSurfaceOpd",
)

const HR_8799_H_MAGNITUDE = 5.280f0
const SPIDERS_ATMOSPHERE_RNG_SEED = UInt64(0x4852_3837_3939)

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

"""Build the shared, provisional atmosphere-to-LLOWFS/SCC viewer graph."""
function graph_definition(
    llowfs_configuration,
    scc_configuration,
    pupil_amplitude,
    atmosphere_step,
)
    pupil_resolution = llowfs_configuration.pupil_resolution
    pupil_resolution == scc_configuration.pupil_resolution ||
        error("LLOWFS and SCC pupil resolutions must match")
    diameter_m = llowfs_configuration.diameter_m
    diameter_m == scc_configuration.diameter_m ||
        error("LLOWFS and SCC telescope diameters must match")
    pupil_opd_schema = llowfs_configuration.pupil_opd_schema
    pupil_opd_schema == scc_configuration.pupil_opd_schema ||
        error("LLOWFS and SCC pupil OPD schemas must match")
    atmosphere = multilayer_atmosphere_opd_node(
        :atmosphere;
        resolution=pupil_resolution,
        telescope_diameter_m=diameter_m,
        central_obstruction_ratio=0.0f0,
        pupil_reflectivity=1.0f0,
        aperture_revision=1,
        r0=0.16f0,
        reference_wavelength_m=500.0f-9,
        L0=25.0f0,
        fractional_cn2=(0.55f0, 0.20f0, 0.15f0, 0.10f0),
        wind_speed=(5.0f0, 7.0f0, 10.0f0, 12.0f0),
        wind_direction_deg=(0.0f0, 75.0f0, 170.0f0, 260.0f0),
        altitude=(0.0f0, 2_000.0f0, 7_000.0f0, 12_000.0f0),
        layer_ids=(:ground, :two_km, :seven_km, :twelve_km),
        atmosphere_step,
        rng_seed=SPIDERS_ATMOSPHERE_RNG_SEED,
        atmosphere_opd_schema=pupil_opd_schema,
        T=Float32,
    )
    return algorithm_graph(
        (
            atmosphere,
            proper_propagation_node(:llowfs_optics, llowfs_configuration),
            proper_propagation_node(:scc_optics, scc_configuration),
        );
        name=:spiders_hr_8799_viewer,
        inputs=(
            graph_input(
                :llowfs_pupil_amplitude,
                :llowfs_optics => :pupil_amplitude,
                pupil_amplitude,
            ),
            graph_input(
                :scc_pupil_amplitude,
                :scc_optics => :pupil_amplitude,
                pupil_amplitude,
            ),
        ),
        outputs=(
            graph_output(
                :atmosphere_opd,
                :atmosphere => :atmosphere_opd,
            ),
            graph_output(
                :llowfs_relative_intensity,
                :llowfs_optics => :output,
            ),
            graph_output(
                :scc_relative_intensity,
                :scc_optics => :output,
            ),
        ),
        links=(
            link(
                :atmosphere => :atmosphere_opd,
                :llowfs_optics => :pupil_opd,
            ),
            link(
                :atmosphere => :atmosphere_opd,
                :scc_optics => :pupil_opd,
            ),
        ),
    )
end

"""
Convert the prescription's unit-amplitude relative intensity to provisional
source-scaled photon-arrival rate. The scale treats each input pupil sample as
one square collecting-area cell and compensates for the final camera
magnification. It does not supply as-built optical throughput or detector
response.
"""
function provisional_photon_rate_scale(profile, configuration, source)
    pupil_sample_pitch_m =
        profile.optics.subaru_pupil_diameter_m /
        configuration.pupil_resolution
    camera_magnification =
        configuration.prescription.parameters.camera_magnification
    source_photon_irradiance_m2_s = photon_irradiance(source)
    return Float32(
        source_photon_irradiance_m2_s * pupil_sample_pitch_m^2 /
        camera_magnification^2,
    )
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
    source = Source(
        band=:H,
        magnitude=HR_8799_H_MAGNITUDE,
        wavelength=profile.wavelength_um * 1.0f-6,
        T=Float32,
    )
    pupil_shape = (
        llowfs_configuration.pupil_resolution,
        llowfs_configuration.pupil_resolution,
    )
    pupil_opd_host = zeros(Float32, pupil_shape)
    pupil_amplitude_host = provisional_spiders_entrance_pupil_amplitude(
        profile,
        pupil_shape[1],
    )
    pupil_amplitude = device_copy(target, pupil_amplitude_host)
    graph = prepare_algorithm_graph(
        graph_definition(
            llowfs_configuration,
            scc_configuration,
            pupil_amplitude,
            inv(Float32(frame_rate)),
        );
        target,
        execution=StreamGraphExecution(),
    )

    llowfs = zeros(Float32, profile.llowfs_output_shape)
    scc_frame = zeros(Float32, profile.scc_output_shape)
    llowfs_photon_rate_scale = provisional_photon_rate_scale(
        profile,
        llowfs_configuration,
        source,
    )
    scc_photon_rate_scale = provisional_photon_rate_scale(
        profile,
        scc_configuration,
        source,
    )
    pair_plan, pair_state, pair_products = prepare_scc_pairing(
        Float32;
        frame_shape=profile.scc_output_shape,
    )
    dm_surface_opd = zeros(Float32, pupil_shape)
    pupil_opd_stream = nothing
    llowfs_stream = nothing
    scc_unfringed_stream = nothing
    scc_difference_stream = nothing
    dm_surface_opd_stream = nothing
    viewer = nothing
    try
        pupil_opd_stream = create_stream(
            STREAM_NAMES[1],
            Float32,
            pupil_shape,
        )
        llowfs_stream = create_stream(
            STREAM_NAMES[2],
            Float32,
            profile.llowfs_output_shape,
        )
        scc_unfringed_stream = create_stream(
            STREAM_NAMES[3],
            Float32,
            profile.scc_output_shape,
        )
        scc_difference_stream = create_stream(
            STREAM_NAMES[4],
            Float32,
            profile.scc_output_shape,
        )
        dm_surface_opd_stream = create_stream(
            STREAM_NAMES[5],
            Float32,
            pupil_shape,
        )
        command = Cmd(String[
            viewer_executable(),
            STREAM_NAMES...,
            "--geometry",
            "2x3",
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

        println(
            "target HR 8799, H magnitude ", HR_8799_H_MAGNITUDE,
            ", photon irradiance ",
            photon_irradiance(source),
            " photons s^-1 m^-2",
        )

        started = time()
        next_frame = started
        frame = 0
        fringed_count = 0
        unfringed_count = 0
        while time() - started < duration && process_running(viewer)
            frame += 1
            step_graph!(graph)
            copyto!(
                pupil_opd_host,
                graph_output(graph, :atmosphere_opd),
            )
            @. pupil_opd_host = ifelse(
                iszero(pupil_amplitude_host),
                0.0f0,
                pupil_opd_host,
            )
            copyto!(
                llowfs,
                graph_output(graph, :llowfs_relative_intensity),
            )
            copyto!(
                scc_frame,
                graph_output(graph, :scc_relative_intensity),
            )
            @. llowfs = max(
                llowfs * llowfs_photon_rate_scale,
                0.0f0,
            )
            @. scc_frame = max(
                scc_frame * scc_photon_rate_scale,
                0.0f0,
            )

            status = spiders_chopper_frame_status(graph, Val(:scc_optics))
            phase = spiders_chopper_phase(status)
            accept_scc_frame!(
                pair_products,
                pair_state,
                pair_plan,
                scc_frame,
                spiders_chopper_sequence(status),
                phase,
            )
            if phase === SpidersFringed
                fringed_count += 1
            else
                unfringed_count += 1
            end

            publish!(pupil_opd_stream, pupil_opd_host)
            publish!(llowfs_stream, llowfs)
            pair_state.have_unfringed &&
                publish!(scc_unfringed_stream, pair_state.unfringed)
            pair_products.valid &&
                publish!(scc_difference_stream, pair_products.difference)
            publish!(dm_surface_opd_stream, dm_surface_opd)

            if frame == 1 || frame % max(1, round(Int, frame_rate)) == 0
                println(
                    "frame ", frame,
                    ", SCC phase ", phase,
                    ", atmosphere OPD range ", extrema(pupil_opd_host),
                    ", LLOWFS range ", extrema(llowfs),
                    ", unfringed SCC range ", extrema(pair_state.unfringed),
                    ", SCC difference range ", extrema(pair_products.difference),
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
        close_stream_noexcept!(dm_surface_opd_stream)
        close_stream_noexcept!(scc_difference_stream)
        close_stream_noexcept!(scc_unfringed_stream)
        close_stream_noexcept!(llowfs_stream)
        close_stream_noexcept!(pupil_opd_stream)
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
