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
        return _Backends.compute_device(AMDGPU.ROCArray(zeros(Float32, 1)))
    elseif name === :cuda
        CUDA.functional() || error("CUDA is not functional on this host")
        CUDA.allowscalar(false)
        return _Backends.compute_device(CUDA.CuArray(zeros(Float32, 1)))
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
            proper_hil=_repository_state(
                joinpath(_ROOT, "..", "AdaptiveOpticsProperHIL.jl"),
            ),
            adaptive_optics_sim=_repository_state(
                joinpath(_ROOT, "..", "AdaptiveOpticsSim.jl"),
            ),
            proper=_repository_state(joinpath(_ROOT, "..", "Proper.jl")),
        ),
    )
end

function _configurations(resolution::Int)
    profile = provisional_spiders_h_regular_profile(Float32)
    llowfs = provisional_spiders_llowfs_configuration(
        profile;
        resolution,
        pupil_resolution=128,
    )
    scc = provisional_spiders_chopped_scc_configuration(
        profile;
        resolution,
        pupil_resolution=128,
    )
    return (; profile, llowfs, scc)
end

function _graph_definition(
    mode::Val{:independent},
    configurations,
    pupil_opd,
    pupil_amplitude,
)
    return algorithm_graph(
        (
            proper_propagation_node(:llowfs, configurations.llowfs),
            proper_propagation_node(:scc, configurations.scc),
        );
        name=:spiders_independent_optics,
        inputs=(
            graph_input(:llowfs_opd, :llowfs => :pupil_opd, pupil_opd),
            graph_input(
                :llowfs_amplitude,
                :llowfs => :pupil_amplitude,
                pupil_amplitude,
            ),
            graph_input(:scc_opd, :scc => :pupil_opd, pupil_opd),
            graph_input(
                :scc_amplitude,
                :scc => :pupil_amplitude,
                pupil_amplitude,
            ),
        ),
        outputs=(
            graph_output(:llowfs, :llowfs => :output),
            graph_output(:scc, :scc => :output),
        ),
    )
end

function _graph_definition(
    mode::Val{:shared},
    configurations,
    pupil_opd,
    pupil_amplitude,
)
    return algorithm_graph(
        (
            provisional_spiders_optical_node(
                :optics,
                configurations.llowfs,
                configurations.scc,
            ),
        );
        name=:spiders_shared_optics,
        inputs=(
            graph_input(:pupil_opd, :optics => :pupil_opd, pupil_opd),
            graph_input(
                :pupil_amplitude,
                :optics => :pupil_amplitude,
                pupil_amplitude,
            ),
        ),
        outputs=(
            graph_output(:llowfs, :optics => :llowfs_output),
            graph_output(:scc, :optics => :scc_output),
        ),
    )
end

function _device_copy(target, source::Matrix{T}) where {T}
    destination = _Backends.allocate_device_array(target, T, size(source)...)
    copyto!(destination, source)
    return destination
end

function _prepare_graph(mode, configurations, target)
    pupil_opd_host = zeros(Float32, 128, 128)
    pupil_amplitude_host = provisional_spiders_entrance_pupil_amplitude(
        configurations.profile,
        128,
    )
    pupil_opd = _device_copy(target, pupil_opd_host)
    pupil_amplitude = _device_copy(target, pupil_amplitude_host)
    graph = prepare_algorithm_graph(
        _graph_definition(
            mode,
            configurations,
            pupil_opd,
            pupil_amplitude,
        );
        target,
        execution=StreamGraphExecution(),
    )
    return (; graph, pupil_opd, pupil_amplitude)
end

function _relative_l2(actual, reference)
    denominator = norm(reference)
    difference = norm(actual .- reference)
    return iszero(denominator) ? difference : difference / denominator
end

function _outputs(graph)
    return (
        llowfs=Array(graph_output(graph, :llowfs)),
        scc=Array(graph_output(graph, :scc)),
    )
end

function _cpu_reference(configurations)
    target = _Backends.HostComputeDevice()
    prepared = _prepare_graph(Val(:shared), configurations, target)
    step_graph!(prepared.graph)
    return _outputs(prepared.graph)
end

function _check_correctness!(prepared, reference, backend::Val)
    reset_graph!(prepared.graph)
    step_graph!(prepared.graph)
    _synchronize(backend)
    actual = _outputs(prepared.graph)
    result = (
        llowfs_relative_l2=_relative_l2(actual.llowfs, reference.llowfs),
        scc_relative_l2=_relative_l2(actual.scc, reference.scc),
        llowfs_max_abs=maximum(abs, actual.llowfs .- reference.llowfs),
        scc_max_abs=maximum(abs, actual.scc .- reference.scc),
    )
    for name in (:llowfs_relative_l2, :scc_relative_l2)
        value = getproperty(result, name)
        isfinite(value) && value <= 1e-4 || error(
            "$name=$value exceeds the 1e-4 correctness limit",
        )
    end
    reset_graph!(prepared.graph)
    return result
end

function _warm!(graph, count::Int)
    for _ in 1:count
        step_graph!(graph)
    end
    return graph
end

@inline function _step!(graph, llowfs_host, scc_host, ::Val{false})
    step_graph!(graph)
    return nothing
end

@inline function _step!(graph, llowfs_host, scc_host, ::Val{true})
    step_graph!(graph)
    copyto!(llowfs_host, graph_output(graph, :llowfs))
    copyto!(scc_host, graph_output(graph, :scc))
    return nothing
end

function _measure_runs!(
    graph,
    llowfs_host,
    scc_host,
    publish::Val,
    samples::Int,
    repetitions::Int,
)
    runs = Vector{Vector{UInt64}}(undef, repetitions)
    for repetition in eachindex(runs)
        latencies = Vector{UInt64}(undef, samples)
        for index in eachindex(latencies)
            started = time_ns()
            _step!(graph, llowfs_host, scc_host, publish)
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
    llowfs = graph_output(graph, :llowfs)
    scc = graph_output(graph, :scc)
    llowfs_host = Matrix{eltype(llowfs)}(undef, size(llowfs))
    scc_host = Matrix{eltype(scc)}(undef, size(scc))
    _step!(graph, llowfs_host, scc_host, publish)
    synchronized_host_bytes = @allocated _step!(
        graph,
        llowfs_host,
        scc_host,
        publish,
    )
    runs = _measure_runs!(
        graph,
        llowfs_host,
        scc_host,
        publish,
        samples,
        repetitions,
    )
    return (; synchronized_host_bytes, summary=_summary(runs), runs_ns=runs)
end

function _benchmark_mode(
    mode::Val{Mode},
    configurations,
    target,
    backend,
    reference,
    samples,
    repetitions,
    warmup,
) where {Mode}
    prepared = _prepare_graph(mode, configurations, target)
    correctness = _check_correctness!(prepared, reference, backend)
    _warm!(prepared.graph, warmup)
    core = _measure_phase!(
        prepared.graph,
        Val(false),
        samples,
        repetitions,
    )
    publish = _measure_phase!(
        prepared.graph,
        Val(true),
        samples,
        repetitions,
    )
    return (; mode=Mode, correctness, phases=(; core, publish))
end

function main(args)
    backend_name = Symbol(length(args) >= 1 ? args[1] : "cpu")
    resolution = length(args) >= 2 ? parse(Int, args[2]) : 512
    samples = length(args) >= 3 ? parse(Int, args[3]) : 100
    repetitions = length(args) >= 4 ? parse(Int, args[4]) : 3
    output_path = length(args) >= 5 ? args[5] : joinpath(
        @__DIR__,
        "results",
        "$(backend_name)_shared_optics_$(resolution).json",
    )
    samples >= 100 || throw(ArgumentError(
        "at least 100 samples per repetition are required for p99",
    ))
    repetitions >= 3 || throw(ArgumentError(
        "at least three independent repetitions are required",
    ))
    resolution >= 512 && iseven(resolution) || throw(ArgumentError(
        "resolution must be even and at least 512",
    ))

    warmup = 10
    BLAS.set_num_threads(1)
    _Backends.set_fft_provider_threads!(1)
    backend = Val(backend_name)
    target = _load_backend(backend_name)
    configurations = _configurations(resolution)
    reference = _cpu_reference(configurations)
    _synchronize(backend)
    workloads = (
        _benchmark_mode(
            Val(:independent),
            configurations,
            target,
            backend,
            reference,
            samples,
            repetitions,
            warmup,
        ),
        _benchmark_mode(
            Val(:shared),
            configurations,
            target,
            backend,
            reference,
            samples,
            repetitions,
            warmup,
        ),
    )
    result = (
        schema_version=1,
        benchmark="spiders_shared_optical_graph",
        boundary=(
            load_model="closed-loop, concurrency one",
            core="one complete synchronized chopped LLOWFS and SCC optical step",
            publish="core plus both complete detector-product copies to host storage",
            excluded="compilation, FFT planning, instantiation, warmup, detector noise, atmosphere evolution, RTC reconstruction, and DM update",
        ),
        parameters=(;
            backend=backend_name,
            resolution,
            pupil_resolution=128,
            samples,
            repetitions,
            warmup,
        ),
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
        println(workload.mode, ':')
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
