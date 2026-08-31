using CUDA
using LinearAlgebra: norm
using Proper
using SpidersProper
using Statistics

CUDA.allowscalar(false)

percentile_ms(values_ns, percentile) =
    round(quantile(values_ns, percentile) / 1e6; digits=3)

to_device(::Nothing) = nothing
to_device(array::AbstractArray) = CuArray(array)

function to_device(entrance::SpidersProper.StaticSpidersEntrance)
    return SpidersProper.StaticSpidersEntrance(to_device(entrance.field))
end

function to_device(entrance::SpidersProper.ResidualSpidersEntrance)
    return SpidersProper.ResidualSpidersEntrance(
        to_device(entrance.field),
        to_device(entrance.residual_nm),
        entrance.sampling_m,
    )
end

function to_device(entrance::SpidersProper.ExternalSpidersEntrance)
    return SpidersProper.ExternalSpidersEntrance(
        to_device(entrance.field),
        to_device(entrance.pupil_amplitude),
        to_device(entrance.opd_m),
        entrance.pupil_sampling_m,
        entrance.opd_sampling_m,
    )
end

function to_device_named_tuple(values::NamedTuple)
    pairs = (name => to_device(getfield(values, name)) for name in keys(values))
    return (; pairs...)
end

function to_device(prepared::PreparedLLOWFS)
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
    device_prepared = PreparedLLOWFS(
        prepared.config,
        prepared.llowfs_config,
        entrance,
        context,
        wavefront,
        to_device_named_tuple(prepared.apertures),
        to_device(prepared.apodizer_internal),
        to_device(prepared.fpm_opd_internal),
        to_device(prepared.lyot_reflection_internal),
        to_device(prepared.resample_scratch),
        to_device(prepared.output),
        prepared.refractive_indices,
    )
    llowfs_propagate!(device_prepared)
    CUDA.synchronize()
    return device_prepared
end

function benchmark_cuda_prepared_llowfs(
    gridsize::Int,
    samples::Int,
    warmup::Int,
    ao_resolution::Int,
    beam_diameter_fraction::Float64,
)
    residual_nm = zeros(Float32, ao_resolution, ao_resolution)
    config = SpidersConfig(
        reference_pinhole=true,
        beam_diameter_fraction=beam_diameter_fraction,
    )
    cpu_prepared = prepare_llowfs(
        1.55e-6,
        gridsize;
        config,
        ao_residual_nm=residual_nm,
        T=Float32,
    )
    cpu_reference = copy(llowfs_propagate!(cpu_prepared))
    prepared = to_device(cpu_prepared)

    for _ in 1:warmup
        llowfs_propagate!(prepared)
    end
    CUDA.synchronize()

    latencies_ns = Vector{UInt64}(undef, samples)
    for index in eachindex(latencies_ns)
        started_ns = time_ns()
        llowfs_propagate!(prepared)
        CUDA.synchronize()
        latencies_ns[index] = time_ns() - started_ns
    end
    sorted_ns = sort(latencies_ns)
    gpu_output = Array(prepared.output)
    difference = gpu_output .- cpu_reference
    relative_l2 = norm(difference) / norm(cpu_reference)

    llowfs_propagate!(prepared)
    CUDA.synchronize()
    host_allocated = @allocated begin
        llowfs_propagate!(prepared)
        CUDA.synchronize()
    end

    println("Prepared SPIDERS LLOWFS PROPER CUDA propagation benchmark")
    println("policy: CPU preparation, device adaptation, compilation, and warmup are outside samples")
    println("boundary: bound device input through synchronized reflected-Lyot sensor-plane intensity")
    println("environment:")
    println("  julia_version: ", VERSION)
    println("  julia_threads: ", Threads.nthreads())
    println("  cuda_runtime: ", CUDA.runtime_version())
    println("  cuda_driver: ", CUDA.driver_version())
    println("  device: ", CUDA.name(CUDA.device()))
    println("  capability: ", CUDA.capability(CUDA.device()))
    println("parameters:")
    println("  wavelength_m: 1.55e-6")
    println("  gridsize: ", gridsize)
    println("  beam_diameter_fraction: ", beam_diameter_fraction)
    println("  ao_resolution: ", ao_resolution)
    println("  precision: Float32")
    println("  warmup: ", warmup)
    println("  samples: ", samples)
    println("steady_state:")
    println("  min_ms: ", round(minimum(sorted_ns) / 1e6; digits=3))
    println("  median_ms: ", percentile_ms(sorted_ns, 0.50))
    samples >= 10 && println("  p90_ms: ", percentile_ms(sorted_ns, 0.90))
    samples >= 20 && println("  p95_ms: ", percentile_ms(sorted_ns, 0.95))
    samples >= 100 && println("  p99_ms: ", percentile_ms(sorted_ns, 0.99))
    println("  max_ms: ", round(maximum(sorted_ns) / 1e6; digits=3))
    println("  mean_frames_per_s: ", round(1e9 / mean(sorted_ns); digits=3))
    println("  host_allocated_bytes_per_call: ", host_allocated)
    println("correctness:")
    println("  intensity_checksum: ", sum(gpu_output))
    println("  max_abs_difference_vs_cpu: ", maximum(abs, difference))
    println("  relative_l2_difference_vs_cpu: ", relative_l2)
    println("  output_sampling_m: ", prop_get_sampling(prepared.wavefront))
    println("raw_latencies_ns: ", join(latencies_ns, ','))
    return nothing
end

function main(args)
    gridsize = length(args) >= 1 ? parse(Int, args[1]) : 512
    samples = length(args) >= 2 ? parse(Int, args[2]) : 100
    warmup = length(args) >= 3 ? parse(Int, args[3]) : 5
    ao_resolution = length(args) >= 4 ? parse(Int, args[4]) : 128
    beam_diameter_fraction = length(args) >= 5 ? parse(Float64, args[5]) : 0.1
    benchmark_cuda_prepared_llowfs(
        gridsize,
        samples,
        warmup,
        ao_resolution,
        beam_diameter_fraction,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
