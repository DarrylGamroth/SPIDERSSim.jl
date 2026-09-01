using CUDA
using Proper
using SpidersProper

include(joinpath(@__DIR__, "benchmark_cuda_graph_llowfs.jl"))

function to_device(prepared::PreparedSpiders)
    entrance = to_device(prepared.entrance)
    field = similar(entrance.field)
    context = Proper.RunContext(typeof(field))
    wavefront = prop_begin!(
        field,
        prepared.config.telescope_diameter_m,
        prop_get_wavelength(prepared.wavefront);
        beam_diam_fraction=prepared.config.beam_diameter_fraction,
        context,
    )
    device_prepared = PreparedSpiders(
        prepared.config,
        entrance,
        context,
        wavefront,
        to_device_named_tuple(prepared.apertures),
        to_device(prepared.apodizer_internal),
        to_device(prepared.fpm_opd_internal),
        to_device(prepared.lyot_internal),
        to_device(prepared.resample_scratch),
        to_device(prepared.output),
    )
    spiders_propagate!(device_prepared)
    CUDA.synchronize()
    return device_prepared
end

function prepare_scc_graph_fixture(
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
        prepare_spiders(
            1.25e-6,
            gridsize;
            config,
            pupil_amplitude,
            pupil_sampling_m=7.92 / 128,
            opd_m,
            T=Float32,
        )
    elseif input_mode === :residual
        residual_nm = zeros(Float32, 128, 128)
        prepare_spiders(
            1.25e-6,
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
        spiders_propagate!(prepared)
    end
    CUDA.synchronize()
    return cpu_prepared, prepared
end

function main_scc(args)
    gridsize = length(args) >= 1 ? parse(Int, args[1]) : 1024
    samples = length(args) >= 2 ? parse(Int, args[2]) : 100
    beam_diameter_fraction = length(args) >= 3 ?
        parse(Float64, args[3]) : 0.1141323837216534
    input_mode = length(args) >= 4 ? Symbol(args[4]) : :external

    cpu_prepared, prepared = prepare_scc_graph_fixture(
        gridsize,
        beam_diameter_fraction,
        input_mode,
    )
    graph = prepare_cuda_graph(spiders_propagate!, prepared; warmup=10)
    correctness = graph_correctness(
        spiders_propagate!,
        cpu_prepared,
        prepared,
        graph,
    )

    for _ in 1:10
        spiders_propagate!(prepared)
    end
    CUDA.synchronize()
    direct_before = measure_direct(spiders_propagate!, prepared, samples)
    graph_latencies = measure_graph(graph, samples)
    direct_after = measure_direct(spiders_propagate!, prepared, samples)

    direct_host_bytes = @allocated begin
        spiders_propagate!(prepared)
        CUDA.synchronize()
    end
    graph_host_bytes = @allocated begin
        launch_cuda_graph!(graph)
        CUDA.synchronize()
    end

    println("Prepared SPIDERS SCC CUDA Graph experiment")
    println("policy: preparation, compilation, capture, instantiation, and warmup are outside samples")
    println("boundary: one synchronized entrance-to-SCC-focus propagation with device-resident inputs and output")
    println("environment:")
    println("  julia_version: ", VERSION)
    println("  cuda_runtime: ", CUDA.runtime_version())
    println("  cuda_driver: ", CUDA.driver_version())
    println("  device: ", CUDA.name(CUDA.device()))
    println("parameters:")
    println("  wavelength_m: 1.25e-6")
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

    GC.@preserve graph CUDA.synchronize()
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_scc(ARGS)
end
