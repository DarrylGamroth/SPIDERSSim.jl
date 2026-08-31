using CUDA
using Profile
using Proper
using SpidersProper

include(joinpath(@__DIR__, "benchmark_cuda_prepared_llowfs.jl"))

function prepare_device_llowfs(
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
    cpu_reference = copy(llowfs_propagate!(cpu_prepared))
    prepared = to_device(cpu_prepared)
    for _ in 1:10
        llowfs_propagate!(prepared)
    end
    CUDA.synchronize()
    return prepared, cpu_reference
end

function measured_host_bytes(prepared)
    llowfs_propagate!(prepared)
    CUDA.synchronize()
    propagation_bytes = @allocated llowfs_propagate!(prepared)
    synchronization_bytes = @allocated CUDA.synchronize()
    combined_bytes = @allocated begin
        llowfs_propagate!(prepared)
        CUDA.synchronize()
    end
    return (; propagation_bytes, synchronization_bytes, combined_bytes)
end

function collect_allocations(prepared)
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=1.0 llowfs_propagate!(prepared)
    CUDA.synchronize()
    return copy(Profile.Allocs.fetch().allocs)
end

function allocation_site(allocation)
    for frame in allocation.stacktrace
        file = string(frame.file)
        if occursin("/CUDA/", file) ||
           occursin("/Proper/", file) ||
           occursin("/SpidersProper/", file) ||
           occursin("spidersProper", file)
            return (string(frame.func), file, frame.line, string(allocation.type))
        end
    end
    for frame in allocation.stacktrace
        file = string(frame.file)
        if frame.line > 0 &&
           !endswith(file, ".c") &&
           !endswith(file, ".h")
            return (string(frame.func), file, frame.line, string(allocation.type))
        end
    end
    return ("unknown", "unknown", 0, string(allocation.type))
end

function print_allocation_summary(allocations; limit::Int=30)
    groups = Dict{Tuple{String,String,Int,String},Tuple{Int,Int}}()
    for allocation in allocations
        site = allocation_site(allocation)
        count, bytes = get(groups, site, (0, 0))
        groups[site] = (count + 1, bytes + allocation.size)
    end
    ranked = sort!(collect(groups); by=entry -> entry.second[2], rev=true)

    println("allocation_profile:")
    println("  sampled_allocations: ", length(allocations))
    println("  sampled_bytes: ", sum(allocation -> allocation.size, allocations; init=0))
    println("  sample_rate: 1.0")
    println("  top_sites:")
    for (site, (count, bytes)) in Iterators.take(ranked, limit)
        function_name, file, line, allocation_type = site
        println(
            "    bytes=", bytes,
            " count=", count,
            " type=", allocation_type,
            " at=", function_name,
            " ", file,
            ":", line,
        )
    end
    return nothing
end

function main(args)
    gridsize = length(args) >= 1 ? parse(Int, args[1]) : 1024
    beam_diameter_fraction = length(args) >= 2 ?
        parse(Float64, args[2]) : 0.1141323837216534
    input_mode = length(args) >= 3 ? Symbol(args[3]) : :external

    prepared, cpu_reference = prepare_device_llowfs(
        gridsize,
        beam_diameter_fraction,
        input_mode,
    )
    println("Prepared SPIDERS LLOWFS CUDA host-allocation profile")
    println("julia_version: ", VERSION)
    println("cuda_runtime: ", CUDA.runtime_version())
    println("device: ", CUDA.name(CUDA.device()))
    println("gridsize: ", gridsize)
    println("beam_diameter_fraction: ", beam_diameter_fraction)
    println("input_mode: ", input_mode)
    println("measured_host_bytes: ", measured_host_bytes(prepared))

    # Warm the allocation-profiler instrumentation before retaining its data.
    collect_allocations(prepared)
    allocations = collect_allocations(prepared)
    print_allocation_summary(allocations)

    gpu_output = Array(prepared.output)
    relative_l2 = sqrt(sum(abs2, gpu_output .- cpu_reference)) /
                  sqrt(sum(abs2, cpu_reference))
    println("correctness:")
    println("  relative_l2_difference_vs_cpu: ", relative_l2)
    println("  output_checksum: ", sum(gpu_output))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
