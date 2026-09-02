using AdaptiveOpticsProperHIL
using AdaptiveOpticsSim
import AMDGPU
import CUDA
using JSON3
using LinearAlgebra
using SPIDERSSim
using Statistics

using AdaptiveOpticsSim.AlgorithmGraphs

const _Backends = AdaptiveOpticsSim.Backends
const _ROOT = normpath(joinpath(@__DIR__, ".."))

function _load_backend(name::Symbol)
    if name === :cpu
        return _Backends.HostComputeDevice()
    elseif name === :amdgpu
        AMDGPU.functional() || error("AMDGPU is not functional on this host")
        AMDGPU.allowscalar(false)
        probe = AMDGPU.ROCArray(zeros(Float32, 1))
        return _Backends.compute_device(probe)
    elseif name === :cuda
        CUDA.functional() || error("CUDA is not functional on this host")
        CUDA.allowscalar(false)
        probe = CUDA.CuArray(zeros(Float32, 1))
        return _Backends.compute_device(probe)
    end
    throw(ArgumentError("backend must be cpu, amdgpu, or cuda"))
end

_synchronize(::Val{:cpu}) = nothing
_synchronize(::Val{:amdgpu}) = AMDGPU.synchronize()
_synchronize(::Val{:cuda}) = CUDA.synchronize()

function _backend_details(::Val{:cpu})
    model = "unknown"
    for line in eachline("/proc/cpuinfo")
        startswith(line, "model name") || continue
        model = strip(last(split(line, ':'; limit=2)))
        break
    end
    return (device=model, runtime="host", driver="kernel")
end

function _backend_details(::Val{:amdgpu})
    return (
        device=string(AMDGPU.device()),
        package=string(pkgversion(AMDGPU)),
        runtime=string(AMDGPU.HIP.runtime_version()),
        fft=string(AMDGPU.rocFFT.version()),
        driver="kernel amdgpu",
    )
end

function _backend_details(::Val{:cuda})
    return (
        device=CUDA.name(CUDA.device()),
        package=string(pkgversion(CUDA)),
        runtime=string(CUDA.runtime_version()),
        driver=string(CUDA.driver_version()),
    )
end

function _repository_state(path::AbstractString)
    revision = readchomp(`git -C $path rev-parse HEAD`)
    status = readchomp(`git -C $path status --short`)
    return (revision, dirty=!isempty(status))
end

function _environment(backend::Val)
    return (
        julia_version=string(VERSION),
        julia_threads=Threads.nthreads(),
        blas_threads=BLAS.get_num_threads(),
        kernel=readchomp(`uname -srvm`),
        active_project=Base.active_project(),
        backend=_backend_details(backend),
        repositories=(
            spiders_sim=_repository_state(_ROOT),
            proper_hil=_repository_state(joinpath(_ROOT, "..", "AdaptiveOpticsProperHIL.jl")),
            adaptive_optics_sim=_repository_state(joinpath(_ROOT, "..", "AdaptiveOpticsSim.jl")),
            proper=_repository_state(joinpath(_ROOT, "..", "Proper.jl")),
        ),
    )
end

function _tilted_opd!(opd::AbstractMatrix{T}, amplitude::T) where {T}
    center = T(size(opd, 2) + 1) / T(2)
    @inbounds for column in axes(opd, 2)
        value = amplitude * (T(column) - center) / T(size(opd, 2))
        for row in axes(opd, 1)
            opd[row, column] = value
        end
    end
    return opd
end

function _graph_definition(configuration, pupil_opd, pupil_amplitude, name)
    return algorithm_graph(
        (proper_propagation_node(:optics, configuration),);
        name,
        inputs=(
            graph_input(:pupil_opd, :optics => :pupil_opd, pupil_opd),
            graph_input(
                :pupil_amplitude,
                :optics => :pupil_amplitude,
                pupil_amplitude,
            ),
        ),
        outputs=(graph_output(:intensity, :optics => :output),),
    )
end

function _device_copy(target, source::Matrix{T}) where {T}
    destination = _Backends.allocate_device_array(target, T, size(source)...)
    copyto!(destination, source)
    return destination
end

function _prepare_graph(
    configuration,
    target,
    execution,
    name::Symbol,
)
    pupil_shape = (configuration.pupil_resolution, configuration.pupil_resolution)
    pupil_opd_host = zeros(Float32, pupil_shape)
    pupil_amplitude_host = ones(Float32, pupil_shape)
    pupil_opd = _device_copy(target, pupil_opd_host)
    pupil_amplitude = _device_copy(target, pupil_amplitude_host)
    graph = prepare_algorithm_graph(
        _graph_definition(configuration, pupil_opd, pupil_amplitude, name);
        target,
        execution,
    )
    return (; graph, pupil_opd, pupil_amplitude)
end

function _cpu_reference(configuration, tilted::Bool)
    prepared = _prepare_graph(
        configuration,
        _Backends.HostComputeDevice(),
        StreamGraphExecution(),
        :cpu_reference,
    )
    if tilted
        _tilted_opd!(prepared.pupil_opd, 2f-7)
    end
    step_graph!(prepared.graph)
    return copy(graph_output(prepared.graph, :intensity))
end

function _relative_l2(actual, reference)
    difference = actual .- reference
    denominator = norm(reference)
    return iszero(denominator) ? norm(difference) : norm(difference) / denominator
end

function _validate_correctness(correctness; relative_l2_limit=1e-4)
    for (name, value) in pairs(correctness)
        isfinite(value) || error("$name correctness result is not finite")
        occursin("relative_l2", String(name)) || continue
        value <= relative_l2_limit || error(
            "$name=$value exceeds relative L2 limit $relative_l2_limit",
        )
    end
    return correctness
end

function _check_correctness!(
    prepared,
    configuration,
    backend::Val,
)
    flat_reference = _cpu_reference(configuration, false)
    tilted_reference = _cpu_reference(configuration, true)

    fill!(prepared.pupil_opd, 0f0)
    _synchronize(backend)
    step_graph!(prepared.graph)
    flat = Array(graph_output(prepared.graph, :intensity))

    tilted_host = zeros(Float32, size(prepared.pupil_opd))
    _tilted_opd!(tilted_host, 2f-7)
    copyto!(prepared.pupil_opd, tilted_host)
    _synchronize(backend)
    step_graph!(prepared.graph)
    tilted = Array(graph_output(prepared.graph, :intensity))

    fill!(prepared.pupil_opd, 0f0)
    _synchronize(backend)
    return _validate_correctness((
        flat_relative_l2=_relative_l2(flat, flat_reference),
        tilted_relative_l2=_relative_l2(tilted, tilted_reference),
        flat_max_abs=maximum(abs, flat .- flat_reference),
        tilted_max_abs=maximum(abs, tilted .- tilted_reference),
    ))
end

function _warm!(graph, count::Int)
    for _ in 1:count
        step_graph!(graph)
    end
    return graph
end

@inline function _step_and_publish!(graph, host_output, ::Val{false})
    step_graph!(graph)
    return nothing
end

@inline function _step_and_publish!(graph, host_output, ::Val{true})
    step_graph!(graph)
    copyto!(host_output, graph_output(graph, :intensity))
    return nothing
end

function _measure_runs!(
    graph,
    host_output,
    ::Val{Publish},
    samples::Int,
    repetitions::Int,
) where {Publish}
    runs = Vector{Vector{UInt64}}(undef, repetitions)
    for repetition in eachindex(runs)
        latencies = Vector{UInt64}(undef, samples)
        for index in eachindex(latencies)
            started = time_ns()
            _step_and_publish!(graph, host_output, Val(Publish))
            latencies[index] = time_ns() - started
        end
        runs[repetition] = latencies
    end
    return runs
end

function _summary(runs::Vector{Vector{UInt64}})
    values = reduce(vcat, runs)
    medians = map(run -> quantile(run, 0.5), runs)
    return (
        samples=length(values),
        repetitions=length(runs),
        minimum_ns=minimum(values),
        p50_ns=quantile(values, 0.50),
        p90_ns=quantile(values, 0.90),
        p95_ns=quantile(values, 0.95),
        p99_ns=quantile(values, 0.99),
        maximum_ns=maximum(values),
        mean_frames_per_s=1e9 / mean(values),
        repetition_p50_min_ns=minimum(medians),
        repetition_p50_max_ns=maximum(medians),
    )
end

function _measure_phase!(graph, publish::Val, samples, repetitions)
    output = graph_output(graph, :intensity)
    host_output = Matrix{eltype(output)}(undef, size(output))
    _step_and_publish!(graph, host_output, publish)
    synchronized_host_bytes = @allocated _step_and_publish!(
        graph,
        host_output,
        publish,
    )
    runs = _measure_runs!(
        graph,
        host_output,
        publish,
        samples,
        repetitions,
    )
    return (; synchronized_host_bytes, summary=_summary(runs), runs_ns=runs)
end

function _configuration(arm::Symbol)
    profile = provisional_spiders_h_regular_profile(Float32)
    if arm === :llowfs
        return provisional_spiders_llowfs_configuration(profile)
    elseif arm === :scc
        return provisional_spiders_scc_configuration(profile)
    end
    throw(ArgumentError("arm must be llowfs or scc"))
end

function _benchmark_arm(
    arm::Symbol,
    target,
    backend::Val{Backend},
    samples::Int,
    repetitions::Int,
    warmup::Int,
) where {Backend}
    configuration = _configuration(arm)
    direct = _prepare_graph(
        configuration,
        target,
        StreamGraphExecution(),
        Symbol(arm, :_direct),
    )
    _warm!(direct.graph, warmup)
    correctness_direct = _check_correctness!(direct, configuration, backend)

    direct_before = _measure_phase!(
        direct.graph,
        Val(false),
        samples,
        repetitions,
    )
    direct_publish = _measure_phase!(
        direct.graph,
        Val(true),
        samples,
        repetitions,
    )

    if Backend === :cpu
        return (
            arm,
            resolution=configuration.resolution,
            pupil_resolution=configuration.pupil_resolution,
            output_shape=(configuration.output_rows, configuration.output_columns),
            correctness=(direct=correctness_direct,),
            phases=(direct=direct_before, direct_publish),
        )
    end

    captured = _prepare_graph(
        configuration,
        target,
        CapturedGraphExecution(),
        Symbol(arm, :_captured),
    )
    _warm!(captured.graph, warmup)
    correctness_captured = _check_correctness!(
        captured,
        configuration,
        backend,
    )
    captured_core = _measure_phase!(
        captured.graph,
        Val(false),
        samples,
        repetitions,
    )
    captured_publish = _measure_phase!(
        captured.graph,
        Val(true),
        samples,
        repetitions,
    )
    direct_after = _measure_phase!(
        direct.graph,
        Val(false),
        samples,
        repetitions,
    )
    return (
        arm,
        resolution=configuration.resolution,
        pupil_resolution=configuration.pupil_resolution,
        output_shape=(configuration.output_rows, configuration.output_columns),
        correctness=(direct=correctness_direct, captured=correctness_captured),
        phases=(
            direct_before,
            captured=captured_core,
            direct_after,
            direct_publish,
            captured_publish,
        ),
    )
end

function main(args)
    backend_name = Symbol(length(args) >= 1 ? args[1] : "cpu")
    samples = length(args) >= 2 ? parse(Int, args[2]) : 100
    repetitions = length(args) >= 3 ? parse(Int, args[3]) : 3
    output_path = length(args) >= 4 ? args[4] : joinpath(
        @__DIR__,
        "results",
        string(backend_name, "_optical_graph.json"),
    )
    samples >= 100 || throw(ArgumentError(
        "at least 100 samples per repetition are required for p99",
    ))
    repetitions >= 3 || throw(ArgumentError(
        "at least three independent repetitions are required",
    ))
    warmup = 10
    BLAS.set_num_threads(1)
    _Backends.set_fft_provider_threads!(1)
    backend = Val(backend_name)
    target = _load_backend(backend_name)
    _synchronize(backend)

    workloads = (
        _benchmark_arm(:llowfs, target, backend, samples, repetitions, warmup),
        _benchmark_arm(:scc, target, backend, samples, repetitions, warmup),
    )
    result = (
        schema_version=1,
        benchmark="spiders_optical_graph",
        boundary=(
            load_model="closed-loop, concurrency one",
            core="one complete synchronized optical graph step",
            publish="core step plus complete output copy to preallocated host storage",
            excluded="compilation, FFT planning, graph capture, instantiation, warmup, detector noise, RTC reconstruction, and DM update",
        ),
        parameters=(; backend=backend_name, samples, repetitions, warmup),
        environment=_environment(backend),
        workloads,
    )
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        JSON3.pretty(io, result)
        write(io, '\n')
    end
    println("wrote ", output_path)
    for workload in workloads
        println(workload.arm, ':')
        for (name, phase) in pairs(workload.phases)
            summary = phase.summary
            println(
                "  ",
                name,
                " p50_ms=", round(summary.p50_ns / 1e6; digits=3),
                " p99_ms=", round(summary.p99_ns / 1e6; digits=3),
                " fps=", round(summary.mean_frames_per_s; digits=2),
                " host_bytes=", phase.synchronized_host_bytes,
            )
        end
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
