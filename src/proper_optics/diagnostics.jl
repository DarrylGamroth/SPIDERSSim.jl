"""
    psf_contrast(result_or_wavefront; reference_peak=nothing)

Compute the radial 25th, 50th, and 75th percentile intensity profiles in
resolution-element (`lambda/D`) coordinates. Supply `reference_peak` to
normalize the profiles as coronagraphic contrast.
"""
function psf_contrast(input; reference_peak::Union{Nothing,Real}=nothing)
    wf = input isa SCCPropagationResult ? input.wavefront : input
    psf = abs2.(centered_field(wf))
    n = prop_get_gridsize(wf)
    axis = (0:(n - 1)) .- n / 2
    radius = sqrt.(reshape(axis, :, 1) .^ 2 .+ reshape(axis, 1, :) .^ 2)
    radii_pixels = 0:(n ÷ 2 - 1)
    percentiles = (0.25, 0.5, 0.75)
    profiles = Matrix{Float64}(undef, length(percentiles), length(radii_pixels))
    for (column, r) in enumerate(radii_pixels)
        annulus = psf[abs.(radius .- r) .<= 0.5]
        for (row, q) in enumerate(percentiles)
            profiles[row, column] = quantile(annulus, q)
        end
    end
    reference_peak === nothing || (profiles ./= reference_peak)
    separation = collect(radii_pixels) .* prop_get_sampling_radians(wf) /
                 prop_get_wavelength(wf) .* wf.beam_diameter_m
    return (separation_lambda_over_d=separation, percentiles=percentiles, intensity=profiles)
end

function poisson_sample(rng::AbstractRNG, lambda::Real)
    lambda >= 0 || throw(DomainError(lambda, "Poisson mean must be nonnegative"))
    iszero(lambda) && return 0
    if lambda < 30
        threshold = exp(-lambda)
        product = 1.0
        count = 0
        while product > threshold
            count += 1
            product *= rand(rng)
        end
        return count - 1
    end
    # The detector model is approximate already; use the standard Gaussian
    # large-count limit to avoid a heavy distribution dependency.
    return max(0, round(Int, lambda + sqrt(lambda) * randn(rng)))
end

"""
    simulate_cred2(result_or_wavefront; kwargs...)

Resample the final intensity onto a C-RED2-like detector, scale it from stellar
magnitude, add photon and read noise, and return integer ADU values.

This preserves the simplified assumptions of the MATLAB helper: fixed band
zero points and throughput, no dark current, saturation, or nonlinearity.
"""
function simulate_cred2(
    input;
    pixel_size_m::Real=15e-6,
    pixels::Integer=512,
    gain_e_per_adu::Real=33e3 / 2^14,
    quantum_efficiency::Real=0.7,
    read_noise_e::Real=25,
    exposure_s::Real=1,
    magnitude::Real=6,
    rng::AbstractRNG=default_rng(),
)
    wf = input isa SCCPropagationResult ? input.wavefront : input
    pixel_size_m > 0 || throw(ArgumentError("pixel_size_m must be positive"))
    pixels > 0 || throw(ArgumentError("pixels must be positive"))
    gain_e_per_adu > 0 || throw(ArgumentError("gain_e_per_adu must be positive"))
    zero_points = (13.5, 5.63, 4.38, 3.22, 2.82, 1.51) .* 1e9
    wavelength = prop_get_wavelength(wf)
    band_index = wavelength < 0.73e-6 ? 1 :
                 wavelength < 0.98e-6 ? 2 :
                 wavelength < 1.1e-6 ? 3 :
                 wavelength < 1.42e-6 ? 4 :
                 wavelength < 1.9e-6 ? 5 : 6
    collecting_area = pi * (wf.beam_diameter_m / 2)^2
    photoelectrons = zero_points[band_index] * collecting_area *
                     10.0^(-0.4 * magnitude) * exposure_s * quantum_efficiency
    optical_throughput = 0.45 * 10 / 17.45
    intensity = abs2.(centered_field(wf)) .* photoelectrons .* optical_throughput
    sampled = prop_magnify(
        intensity,
        prop_get_sampling(wf) / pixel_size_m,
        pixels;
        CONSERVE=true,
    )
    output = Matrix{Int}(undef, size(sampled))
    for index in eachindex(sampled)
        electrons = poisson_sample(rng, max(0, sampled[index])) + read_noise_e * randn(rng)
        output[index] = floor(Int, electrons / gain_e_per_adu)
    end
    return output
end

"""
    cdi_otf(result_or_wavefront; pupil_diameter=0.9,
            reference_separation=1.575, position_angle_rad=pi/2)

Compute the SCC optical-transfer-function decomposition and coherent
differential-imaging estimates used by the MATLAB model.
"""
function cdi_otf(
    input;
    pupil_diameter::Real=0.9,
    reference_separation::Real=1.575,
    position_angle_rad::Real=pi / 2,
)
    wf = input isa SCCPropagationResult ? input.wavefront : input
    n = prop_get_gridsize(wf)
    sampling = prop_get_sampling(wf)
    angular_sampling = prop_get_sampling_radians(wf)
    wavelength = prop_get_wavelength(wf)
    pinhole_diameter = pupil_diameter / 0.9 / 43 * 24 / 14
    pupil_sampling = wavelength / wf.beam_diameter_m / angular_sampling / n
    xc = reference_separation * cos(position_angle_rad) / pupil_sampling
    yc = reference_separation * sin(position_angle_rad) / pupil_sampling

    mask_pinhole = prop_ellipse(
        wf,
        pinhole_diameter / 2 / pupil_sampling * sampling,
        pinhole_diameter / 2 / pupil_sampling * sampling,
    )
    mask_pinhole_background = prop_ellipse(
        wf,
        2 * pinhole_diameter / pupil_sampling * sampling,
        2 * pinhole_diameter / pupil_sampling * sampling,
    ) .- mask_pinhole
    mask_center = prop_ellipse(
        wf,
        pupil_diameter / pupil_sampling * sampling,
        pupil_diameter / pupil_sampling * sampling,
    )
    mask_side = prop_ellipse(
        wf,
        pupil_diameter / 2 / pupil_sampling * sampling,
        pupil_diameter / 2 / pupil_sampling * sampling,
        xc * sampling,
        yc * sampling,
    )

    intensity = abs2.(centered_field(wf))
    otf = fftshift(fft(intensity))
    i0 = abs.(ifft(ifftshift(otf .* mask_center)))
    side_lobe = circshift(otf .* mask_side, round.(Int, (-yc, -xc)))
    i_minus = abs.(ifft(ifftshift(side_lobe)))
    visibility_map = zeros(Float64, size(i0))
    nonzero = i0 .> 0
    visibility_map[nonzero] .= 2 .* i_minus[nonzero] ./ i0[nonzero]

    background_plane = otf .* mask_pinhole_background
    # MATLAB defines an ordering for complex median values; Julia does not.
    # The original code includes the zero-valued pixels outside the annulus, so
    # the real and imaginary component medians reproduce its zero-dominated
    # background estimate without introducing a complex ordering convention.
    pinhole_background = complex(median(real.(background_plane)), median(imag.(background_plane)))
    reference_intensity = abs.(ifft(ifftshift((otf .- pinhole_background) .* mask_pinhole)))
    planet = zeros(Float64, size(i0))
    valid_reference = reference_intensity .> 0
    planet[valid_reference] .= i0[valid_reference] .-
                              i_minus[valid_reference] .^ 2 ./ reference_intensity[valid_reference] .-
                              reference_intensity[valid_reference]
    planet .= max.(planet, 0)
    planet_visibility = max.(i0 .* (1 .- visibility_map), 0)
    visibility = 2 * sum(i_minus) / sum(i0)
    extinction = i0 ./ planet
    extinction_from_visibility = i0 ./ planet_visibility
    return (
        otf=otf,
        mtf=abs.(otf),
        reference_intensity=reference_intensity,
        modulated_intensity=i_minus,
        unmodulated_intensity=i0,
        planet_intensity=planet,
        planet_intensity_from_visibility=planet_visibility,
        visibility_map=visibility_map,
        visibility=visibility,
        extinction=extinction,
        extinction_from_visibility=extinction_from_visibility,
    )
end
