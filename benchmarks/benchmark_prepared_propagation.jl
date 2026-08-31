using FFTW
using Proper
using SpidersProper
using Statistics

percentile_ms(values_ns, percentile) =
    round(quantile(values_ns, percentile) / 1e6; digits=3)

function allocation_bytes(prepared)
    spiders_propagate!(prepared)
    return @allocated spiders_propagate!(prepared)
end

function benchmark_prepared(
    gridsize::Int,
    samples::Int,
    warmup::Int,
    ao_resolution::Int,
)
    residual_nm = zeros(Float32, ao_resolution, ao_resolution)
    config = SpidersConfig(reference_pinhole=true)
    preparation = @timed prepare_spiders(
        1.25e-6,
        gridsize;
        config,
        ao_residual_nm=residual_nm,
    )
    prepared = preparation.value

    for _ in 1:warmup
        spiders_propagate!(prepared)
    end
    latencies_ns = Vector{UInt64}(undef, samples)
    for index in eachindex(latencies_ns)
        started_ns = time_ns()
        spiders_propagate!(prepared)
        latencies_ns[index] = time_ns() - started_ns
    end
    sorted_ns = sort(latencies_ns)

    println("Prepared SPIDERS PROPER propagation benchmark")
    println("policy: preparation and compilation are outside warmed closed-loop samples")
    println("boundary: update-bound-input to centered intensity in a caller-stable output buffer")
    println("environment:")
    println("  julia_version: ", VERSION)
    println("  julia_threads: ", Threads.nthreads())
    println("  fftw_threads: ", FFTW.get_num_threads())
    println("  cpu: ", isempty(Sys.cpu_info()) ? "unknown" : Sys.cpu_info()[1].model)
    println("  proper_source: ", pathof(Proper))
    println("parameters:")
    println("  gridsize: ", gridsize)
    println("  ao_resolution: ", ao_resolution)
    println("  warmup: ", warmup)
    println("  samples: ", samples)
    println("preparation:")
    println("  elapsed_s: ", round(preparation.time; digits=6))
    println("  compile_s: ", round(preparation.compile_time; digits=6))
    println("  allocated_bytes: ", preparation.bytes)
    println("steady_state:")
    println("  min_ms: ", round(minimum(sorted_ns) / 1e6; digits=3))
    println("  median_ms: ", percentile_ms(sorted_ns, 0.50))
    samples >= 10 && println("  p90_ms: ", percentile_ms(sorted_ns, 0.90))
    samples >= 20 && println("  p95_ms: ", percentile_ms(sorted_ns, 0.95))
    samples >= 100 && println("  p99_ms: ", percentile_ms(sorted_ns, 0.99))
    println("  max_ms: ", round(maximum(sorted_ns) / 1e6; digits=3))
    println("  mean_frames_per_s: ", round(1e9 / mean(sorted_ns); digits=3))
    println("  allocated_bytes_per_call: ", allocation_bytes(prepared))
    println("correctness:")
    println("  intensity_checksum: ", sum(spiders_intensity(prepared)))
    println("  output_sampling_m: ", prop_get_sampling(prepared.wavefront))
    println("raw_latencies_ns: ", join(latencies_ns, ','))
    return nothing
end

function main(args)
    gridsize = length(args) >= 1 ? parse(Int, args[1]) : 512
    samples = length(args) >= 2 ? parse(Int, args[2]) : 20
    warmup = length(args) >= 3 ? parse(Int, args[3]) : 3
    ao_resolution = length(args) >= 4 ? parse(Int, args[4]) : 128
    benchmark_prepared(gridsize, samples, warmup, ao_resolution)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
