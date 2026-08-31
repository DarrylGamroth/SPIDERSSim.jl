using FFTW
using Proper
using SpidersProper
using Statistics

function percentile_ms(values_ns::Vector{UInt64}, percentile::Float64)
    return round(quantile(values_ns, percentile) / 1e6; digits=3)
end

function benchmark_propagation(
    gridsize::Int,
    samples::Int,
    warmup::Int,
    ao_resolution::Int,
)
    gridsize > 0 || throw(ArgumentError("gridsize must be positive"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    warmup >= 0 || throw(ArgumentError("warmup must be non-negative"))
    ao_resolution > 0 || throw(ArgumentError("ao_resolution must be positive"))

    wavelength_m = 1.25e-6
    ao_residual_nm = zeros(Float32, ao_resolution, ao_resolution)
    config = SpidersConfig(reference_pinhole=true)

    # This deliberately times the public, unprepared API. Static FITS reads,
    # mask construction, propagation workspaces, and the returned wavefront
    # are all inside the operation boundary.
    first_use = @timed spiders_proper(
        wavelength_m,
        gridsize;
        config,
        ao_residual_nm,
    )

    result = first_use.value
    for _ in 1:warmup
        result = spiders_proper(
            wavelength_m,
            gridsize;
            config,
            ao_residual_nm,
        )
    end

    latencies_ns = Vector{UInt64}(undef, samples)
    for sample_index in eachindex(latencies_ns)
        started_ns = time_ns()
        result = spiders_proper(
            wavelength_m,
            gridsize;
            config,
            ao_residual_nm,
        )
        latencies_ns[sample_index] = time_ns() - started_ns
    end

    allocated_bytes = @allocated result = spiders_proper(
        wavelength_m,
        gridsize;
        config,
        ao_residual_nm,
    )
    intensity_checksum = sum(spiders_intensity(result))
    sorted_ns = sort(latencies_ns)

    println("SPIDERS PROPER propagation benchmark")
    println("policy: first propagation after import is separate from warmed closed-loop samples")
    println("boundary: spiders_proper call, including FITS reads, asset/workspace construction, and returned wavefront")
    println("environment:")
    println("  julia_version: ", VERSION)
    println("  julia_threads: ", Threads.nthreads())
    println("  fftw_threads: ", FFTW.get_num_threads())
    println("  cpu: ", isempty(Sys.cpu_info()) ? "unknown" : Sys.cpu_info()[1].model)
    println("  kernel: ", Sys.KERNEL)
    println("  proper_version: ", Base.pkgversion(Proper))
    println("  proper_source: ", pathof(Proper))
    println("parameters:")
    println("  wavelength_m: ", wavelength_m)
    println("  gridsize: ", gridsize)
    println("  ao_resolution: ", ao_resolution)
    println("  ao_input: deterministic zero residual, Float32 nanometers")
    println("  warmup: ", warmup)
    println("  samples: ", samples)
    println("first_use:")
    println("  elapsed_s: ", round(first_use.time; digits=6))
    println("  compile_s: ", round(first_use.compile_time; digits=6))
    println("  gc_s: ", round(first_use.gctime; digits=6))
    println("  allocated_bytes: ", first_use.bytes)
    println("steady_state:")
    println("  min_ms: ", round(minimum(sorted_ns) / 1e6; digits=3))
    println("  median_ms: ", percentile_ms(sorted_ns, 0.50))
    samples >= 10 && println("  p90_ms: ", percentile_ms(sorted_ns, 0.90))
    samples >= 20 && println("  p95_ms: ", percentile_ms(sorted_ns, 0.95))
    samples >= 100 && println("  p99_ms: ", percentile_ms(sorted_ns, 0.99))
    println("  max_ms: ", round(maximum(sorted_ns) / 1e6; digits=3))
    println("  mean_frames_per_s: ", round(1e9 / mean(sorted_ns); digits=3))
    println("  allocated_bytes_per_call_probe: ", allocated_bytes)
    println("correctness:")
    println("  intensity_checksum: ", intensity_checksum)
    println("  output_sampling_m: ", prop_get_sampling(result.wavefront))
    println("raw_latencies_ns: ", join(latencies_ns, ','))
    return nothing
end

function main(args)
    gridsize = length(args) >= 1 ? parse(Int, args[1]) : 512
    samples = length(args) >= 2 ? parse(Int, args[2]) : 20
    warmup = length(args) >= 3 ? parse(Int, args[3]) : 3
    ao_resolution = length(args) >= 4 ? parse(Int, args[4]) : 128
    benchmark_propagation(gridsize, samples, warmup, ao_resolution)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
