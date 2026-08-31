using FFTW
using Proper
using SpidersProper
using Statistics

percentile_ms(values_ns, percentile) =
    round(quantile(values_ns, percentile) / 1e6; digits=3)

function allocation_bytes(prepared)
    llowfs_propagate!(prepared)
    return @allocated llowfs_propagate!(prepared)
end

function benchmark_prepared_llowfs(
    gridsize::Int,
    samples::Int,
    warmup::Int,
    ao_resolution::Int,
    ::Type{T},
    beam_diameter_fraction::Float64,
) where {T<:AbstractFloat}
    residual_nm = zeros(T, ao_resolution, ao_resolution)
    config = SpidersConfig(
        reference_pinhole=true,
        beam_diameter_fraction=beam_diameter_fraction,
    )
    preparation = @timed prepare_llowfs(
        1.55e-6,
        gridsize;
        config,
        ao_residual_nm=residual_nm,
        T,
    )
    prepared = preparation.value

    for _ in 1:warmup
        llowfs_propagate!(prepared)
    end
    latencies_ns = Vector{UInt64}(undef, samples)
    for index in eachindex(latencies_ns)
        started_ns = time_ns()
        llowfs_propagate!(prepared)
        latencies_ns[index] = time_ns() - started_ns
    end
    sorted_ns = sort(latencies_ns)

    println("Prepared SPIDERS LLOWFS PROPER propagation benchmark")
    println("policy: preparation and compilation are outside warmed closed-loop samples")
    println("boundary: bound entrance input through reflected Lyot relay to centered sensor-plane intensity")
    println("environment:")
    println("  julia_version: ", VERSION)
    println("  julia_threads: ", Threads.nthreads())
    println("  fftw_threads: ", FFTW.get_num_threads())
    println("  cpu: ", isempty(Sys.cpu_info()) ? "unknown" : Sys.cpu_info()[1].model)
    println("  proper_source: ", pathof(Proper))
    println("parameters:")
    println("  wavelength_m: 1.55e-6")
    println("  gridsize: ", gridsize)
    println("  beam_diameter_fraction: ", beam_diameter_fraction)
    println("  ao_resolution: ", ao_resolution)
    println("  precision: ", T)
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
    println("  intensity_checksum: ", sum(llowfs_intensity(prepared)))
    println("  output_sampling_m: ", prop_get_sampling(prepared.wavefront))
    println("raw_latencies_ns: ", join(latencies_ns, ','))
    return nothing
end

function main(args)
    gridsize = length(args) >= 1 ? parse(Int, args[1]) : 512
    samples = length(args) >= 2 ? parse(Int, args[2]) : 20
    warmup = length(args) >= 3 ? parse(Int, args[3]) : 3
    ao_resolution = length(args) >= 4 ? parse(Int, args[4]) : 128
    T = length(args) >= 5 ? getproperty(Base, Symbol(args[5])) : Float32
    beam_diameter_fraction = length(args) >= 6 ? parse(Float64, args[6]) : 0.1
    T in (Float32, Float64) || throw(ArgumentError(
        "precision must be Float32 or Float64"))
    benchmark_prepared_llowfs(
        gridsize,
        samples,
        warmup,
        ao_resolution,
        T,
        beam_diameter_fraction,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
