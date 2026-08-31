using CUDA
using LinearAlgebra: norm
using Proper
using SpidersProper
using Statistics

include(joinpath(@__DIR__, "benchmark_cuda_prepared_llowfs.jl"))

function prepare_graph_fixture(
    gridsize::Int,
    beam_diameter_fraction::Float64,
    input_mode::Symbol,
)
    config = SpidersConfig(
        reference_pinhole=true,
        beam_diameter_fraction=beam_diameter_fraction,
    )
    cpu_prepared = if input_mode === :external
        pupil_amplitude = ones(Float32, 128, 128)
        opd_m = zeros(Float32, 128, 128)
        prepare_llowfs(
            1.55e-6,
            gridsize;
            config,
            pupil_amplitude,
            pupil_sampling_m=7.92 / 128,
            opd_m,
            T=Float32,
        )
    elseif input_mode === :residual
        residual_nm = zeros(Float32, 128, 128)
        prepare_llowfs(
            1.55e-6,
            gridsize;
            config,
            ao_residual_nm=residual_nm,
            T=Float32,
        )
    else
        throw(ArgumentError("input mode must be external or residual"))
    end
    prepared = to_device(cpu_prepared)
    for _ in 1:10
        llowfs_propagate!(prepared)
    end
    CUDA.synchronize()
    return cpu_prepared, prepared
end

function capture_propagation(prepared)
    graph = CUDA.capture() do
        llowfs_propagate!(prepared)
    end
    executable = CUDA.instantiate(graph)
    for _ in 1:10
        CUDA.launch(executable)
    end
    CUDA.synchronize()
    return graph, executable
end

function measure_direct(prepared, samples::Int)
    latencies_ns = Vector{UInt64}(undef, samples)
    for index in eachindex(latencies_ns)
        started_ns = time_ns()
        llowfs_propagate!(prepared)
        CUDA.synchronize()
        latencies_ns[index] = time_ns() - started_ns
    end
    return latencies_ns
end

function measure_graph(executable, samples::Int)
    latencies_ns = Vector{UInt64}(undef, samples)
    for index in eachindex(latencies_ns)
        started_ns = time_ns()
        CUDA.launch(executable)
        CUDA.synchronize()
        latencies_ns[index] = time_ns() - started_ns
    end
    return latencies_ns
end

function summarize_latencies(name, values_ns)
    println(name, ':')
    println("  samples: ", length(values_ns))
    println("  min_ms: ", round(minimum(values_ns) / 1e6; digits=3))
    println("  median_ms: ", percentile_ms(values_ns, 0.50))
    println("  p90_ms: ", percentile_ms(values_ns, 0.90))
    println("  p95_ms: ", percentile_ms(values_ns, 0.95))
    length(values_ns) >= 100 &&
        println("  p99_ms: ", percentile_ms(values_ns, 0.99))
    println("  max_ms: ", round(maximum(values_ns) / 1e6; digits=3))
    println("  mean_frames_per_s: ", round(1e9 / mean(values_ns); digits=3))
    return nothing
end

function relative_l2(device_output, cpu_output)
    gpu_output = Array(device_output)
    return norm(gpu_output .- cpu_output) / norm(cpu_output)
end

function graph_correctness(cpu_prepared, prepared, executable)
    cpu_input, device_input, probe = if hasproperty(
        cpu_prepared.entrance,
        :opd_m,
    )
        (cpu_prepared.entrance.opd_m, prepared.entrance.opd_m, 10f-9)
    else
        (
            cpu_prepared.entrance.residual_nm,
            prepared.entrance.residual_nm,
            10f0,
        )
    end

    fill!(cpu_input, 0f0)
    flat_reference = copy(llowfs_propagate!(cpu_prepared))
    fill!(device_input, 0f0)
    CUDA.launch(executable)
    CUDA.synchronize()
    flat_relative_l2 = relative_l2(prepared.output, flat_reference)

    # Change values in the existing input buffers without changing their pointers.
    # A reusable graph must observe this update on every launch.
    fill!(cpu_input, probe)
    probe_reference = copy(llowfs_propagate!(cpu_prepared))
    fill!(device_input, probe)
    CUDA.launch(executable)
    CUDA.synchronize()
    probe_relative_l2 = relative_l2(prepared.output, probe_reference)

    fill!(cpu_input, 0f0)
    fill!(device_input, 0f0)
    return (; flat_relative_l2, probe_relative_l2)
end

function main(args)
    gridsize = length(args) >= 1 ? parse(Int, args[1]) : 1024
    samples = length(args) >= 2 ? parse(Int, args[2]) : 100
    beam_diameter_fraction = length(args) >= 3 ?
        parse(Float64, args[3]) : 0.1141323837216534
    input_mode = length(args) >= 4 ? Symbol(args[4]) : :external

    cpu_prepared, prepared = prepare_graph_fixture(
        gridsize,
        beam_diameter_fraction,
        input_mode,
    )
    graph, executable = capture_propagation(prepared)
    correctness = graph_correctness(cpu_prepared, prepared, executable)

    # Rewarm after the probe and bracket the graph samples with ordinary calls.
    for _ in 1:10
        llowfs_propagate!(prepared)
    end
    CUDA.synchronize()
    direct_before = measure_direct(prepared, samples)
    graph_latencies = measure_graph(executable, samples)
    direct_after = measure_direct(prepared, samples)

    direct_host_bytes = @allocated begin
        llowfs_propagate!(prepared)
        CUDA.synchronize()
    end
    graph_host_bytes = @allocated begin
        CUDA.launch(executable)
        CUDA.synchronize()
    end

    println("Prepared SPIDERS LLOWFS CUDA Graph experiment")
    println("policy: preparation, compilation, capture, instantiation, and warmup are outside samples")
    println("boundary: one synchronized propagation with device-resident inputs and output")
    println("environment:")
    println("  julia_version: ", VERSION)
    println("  cuda_runtime: ", CUDA.runtime_version())
    println("  cuda_driver: ", CUDA.driver_version())
    println("  device: ", CUDA.name(CUDA.device()))
    println("parameters:")
    println("  gridsize: ", gridsize)
    println("  beam_diameter_fraction: ", beam_diameter_fraction)
    println("  input_mode: ", input_mode)
    println("  samples_per_phase: ", samples)
    println("correctness:")
    println("  flat_relative_l2_vs_cpu: ", correctness.flat_relative_l2)
    println("  changed_input_relative_l2_vs_cpu: ", correctness.probe_relative_l2)
    println("allocations:")
    println("  direct_host_bytes_per_call: ", direct_host_bytes)
    println("  graph_host_bytes_per_call: ", graph_host_bytes)
    summarize_latencies("direct_before", direct_before)
    summarize_latencies("graph", graph_latencies)
    summarize_latencies("direct_after", direct_after)
    println("raw_direct_before_ns: ", join(direct_before, ','))
    println("raw_graph_ns: ", join(graph_latencies, ','))
    println("raw_direct_after_ns: ", join(direct_after, ','))

    # Keep graph objects alive until all launches and reporting are complete.
    GC.@preserve graph executable CUDA.synchronize()
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
