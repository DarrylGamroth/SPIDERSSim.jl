"""
    propagate_scc_movie(wavelength_m, ao_residuals_nm; gridsize=1024,
                  configuration=SCCPropagationConfiguration(), final_sampling_m=11.513e-6)

Propagate a cube of AO residual frames through SPIDERS and return a cube of SCC
intensities. The AO cube is indexed `(y, x, time)` and expressed in nanometers.

Unlike the MATLAB routine, the number of frames comes from the input cube and
is not fixed at 1,000. Set `final_sampling_m=nothing` to retain PROPER's native
final sampling and output grid size.
"""
function propagate_scc_movie(
    wavelength_m::Real,
    ao_residuals_nm::AbstractArray{<:Real,3};
    gridsize::Integer=1024,
    configuration::SCCPropagationConfiguration=SCCPropagationConfiguration(reference_pinhole=true),
    ao_sampling_m::Real=0.045 / 7.73 * configuration.telescope_diameter_m,
    final_sampling_m::Union{Nothing,Real}=11.513e-6,
    output_size::Integer=gridsize,
    rng::AbstractRNG=default_rng(),
)
    frames = size(ao_residuals_nm, 3)
    frames > 0 || throw(ArgumentError("ao_residuals_nm must contain at least one frame"))
    output_size > 0 || throw(ArgumentError("output_size must be positive"))
    output = Array{Float64}(undef, output_size, output_size, frames)
    for frame_index in 1:frames
        result = propagate_scc(
            wavelength_m,
            gridsize;
            configuration=configuration,
            ao_residual_nm=@view(ao_residuals_nm[:, :, frame_index]),
            ao_sampling_m=ao_sampling_m,
            rng=rng,
        )
        field = scc_field(result)
        if final_sampling_m === nothing
            size(field) == (output_size, output_size) || throw(ArgumentError(
                "output_size must equal gridsize when final_sampling_m=nothing",
            ))
            @views output[:, :, frame_index] .= abs2.(field)
        else
            magnification = prop_get_sampling(result.wavefront) / final_sampling_m
            resampled = prop_magnify(field, magnification, output_size; CONSERVE=true)
            @views output[:, :, frame_index] .= abs2.(resampled)
        end
    end
    return output
end
