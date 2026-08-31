"""
    subaru_pupil!(wf; obscuration_ratio=0.29, include_spiders=true)

Apply an analytic Subaru pupil, including its central obscuration and optional
secondary-mirror supports, to `wf`. The entrance field must already have been
normalized if unobscured-pupil normalization is required.
"""
function subaru_pupil!(
    wf::Proper.WaveFront;
    obscuration_ratio::Real=0.29,
    include_spiders::Bool=true,
)
    pupil_radius = wf.beam_diameter_m / 2
    prop_circular_obscuration(wf, obscuration_ratio * pupil_radius)
    include_spiders || return wf

    for polygon in subaru_spider_polygons(pupil_radius; margin=1)
        mask = prop_irregular_polygon(wf, polygon.x, polygon.y; DARK=true)
        prop_multiply(wf, mask)
    end
    return wf
end

function subaru_spider_polygons(pupil_radius::Real; margin::Real)
    slope = tand(51.5)
    b2 = -30.7 / 273 * 2 * pupil_radius
    b1 = -22.3 / 273 * 2 * pupil_radius
    thickness = abs(b1 - b2)
    xpositive = pupil_radius .* [0.0, 0.0, 1.0, 1.0]
    ypositive = [b2, b1, slope * pupil_radius + b1, slope * pupil_radius + b2]
    ypositive .+= 0.5 * thickness * (margin - 1) .* [-1, 1, 1, -1]
    ynegative = -ypositive
    xnegative = -xpositive
    return (
        (x=xpositive, y=ypositive),
        (x=xpositive, y=ynegative),
        (x=xnegative, y=-ypositive),
        (x=xnegative, y=-ynegative),
    )
end

"""
    radial_apodizer(coefficients, gridsize, pupil_radius, obscuration_ratio, sampling_m)

Return the normalized centro-symmetric polynomial apodizer used by the MATLAB
SPIDERS model. Coefficients are ordered from constant term upward.
"""
function radial_apodizer(
    coefficients::AbstractVector{<:Real},
    gridsize::Integer,
    pupil_radius::Real,
    obscuration_ratio::Real,
    sampling_m::Real,
)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    pupil_radius > 0 || throw(ArgumentError("pupil_radius must be positive"))
    isempty(coefficients) && throw(ArgumentError("coefficients must not be empty"))

    axis = ((0:(n - 1)) .- n / 2) .* sampling_m
    radius = sqrt.(reshape(axis, :, 1) .^ 2 .+ reshape(axis, 1, :) .^ 2) ./ pupil_radius
    profile = zeros(promote_type(Float64, eltype(coefficients), typeof(sampling_m)), n, n)
    for coefficient in Iterators.reverse(coefficients)
        @. profile = profile * radius + coefficient
    end
    annulus = (radius .>= obscuration_ratio) .& (radius .<= 1)
    any(annulus) || throw(ArgumentError("the sampled grid contains no apodizer pupil annulus"))
    scale = maximum(profile[annulus])
    isfinite(scale) && !iszero(scale) || throw(DomainError(scale, "apodizer normalization is not finite and nonzero"))
    profile ./= scale
    profile .= max.(profile, zero(eltype(profile)))
    return profile
end

"""
    tilt_gaussian_fpm(wf; radius_m, tilt_rad, gaussian_phase=7.8,
                      charge=0, piston_rad=2pi, levels=16,
                      reference_wavelength_m=1.59e-6, center_pupils=true)

Return the centered focal-plane phase map, in radians, for the quantized
Tilt-Gaussian or Tilt-Gaussian-Vortex SPIDERS mask.
"""
function tilt_gaussian_fpm(
    wf::Proper.WaveFront;
    radius_m::Real,
    tilt_rad::Real,
    gaussian_phase::Real=7.8,
    parabola_x::Real=0,
    parabola_y::Real=0,
    charge::Real=0,
    piston_rad::Real=2pi,
    pyramid_apex_m::Real=0,
    levels::Integer=16,
    reference_wavelength_m::Real=1.59e-6,
    center_pupils::Bool=true,
)
    radius_m > 0 || throw(ArgumentError("radius_m must be positive"))
    levels >= 0 || throw(ArgumentError("levels must be nonnegative"))

    n = prop_get_gridsize(wf)
    sampling = prop_get_sampling(wf)
    wavelength = prop_get_wavelength(wf)
    angular_sampling = prop_get_sampling_radians(wf)
    resolution_element_m = wavelength / wf.beam_diameter_m / angular_sampling * sampling
    axis = ((0:(n - 1)) .- n / 2) .* sampling
    xx = reshape(axis, 1, :)
    yy = reshape(axis, :, 1)
    radius = sqrt.(xx .^ 2 .+ yy .^ 2)
    theta = atan.(-xx, -yy)

    core_mask = prop_ellipse(wf, radius_m, radius_m)
    vortex_mask = prop_ellipse(wf, radius_m, radius_m; DARK=true)
    tilt = (-2pi * tilt_rad / wavelength) .* yy

    pyramid = zeros(Float64, n, n)
    if !iszero(pyramid_apex_m)
        p1 = (2pi * pyramid_apex_m / wavelength) .* radius .* cos.(theta)
        p2 = (2pi * pyramid_apex_m / wavelength) .* radius .* cos.(theta .- pi / 4)
        p3 = (2pi * pyramid_apex_m / wavelength) .* radius .* cos.(theta .+ pi / 4)
        first = (-pi / 3 .<= theta) .& (theta .< pi / 3)
        second = (pi / 3 .<= theta) .& (theta .< pi)
        third = (-pi .<= theta) .& (theta .< -pi / 3)
        pyramid[first] .= p1[first]
        pyramid[second] .= p2[second]
        pyramid[third] .= p3[third]
    end

    gaussian_sigma = 2 * resolution_element_m
    gaussian = gaussian_phase .* exp.(-0.5 .* (radius ./ gaussian_sigma) .^ 2)
    parabola = @. -parabola_x / radius_m^2 * xx^2 - parabola_y / radius_m^2 * yy^2
    vortex = charge .* theta .+ piston_rad
    core = tilt .+ parabola .+ gaussian .+ pyramid
    illuminated_core = core_mask .> 0
    core .-= minimum(core[illuminated_core])
    phase = vortex_mask .* vortex .+ core_mask .* core

    if levels > 0
        height_step = reference_wavelength_m / wavelength * 2pi / levels
        phase .= height_step .* floor.(mod.(phase, 2pi) ./ height_step)
    end
    center_pupils && (phase .-= tilt ./ 2)
    return phase
end

"""
    LyotStop(amplitude_transmission, amplitude_reflection,
             main_pupil_transmission, pinhole_transmission)

Contain the centered amplitude coefficients of the SPIDERS reflective Lyot
stop and its incident-flux transmission metrics. For the ideal lossless stop,
the coefficients obey `abs2(amplitude_transmission) +
abs2(amplitude_reflection) == 1` at every sample.
"""
struct LyotStop{A<:AbstractMatrix,T<:Real}
    amplitude_transmission::A
    amplitude_reflection::A
    main_pupil_transmission::T
    pinhole_transmission::T
end

"""Return the centered complex field transmitted by a [`LyotStop`](@ref)."""
transmitted_lyot_field(wf::Proper.WaveFront, stop::LyotStop) =
    centered_field(wf) .* stop.amplitude_transmission

"""Return the centered complex field reflected by a [`LyotStop`](@ref)."""
reflected_lyot_field(wf::Proper.WaveFront, stop::LyotStop) =
    centered_field(wf) .* stop.amplitude_reflection

"""
    scc_lyot_stop(wf, pupil_diameter_m; kwargs...)

Construct the SPIDERS SCC Lyot stop without applying it to `wf`.
"""
function scc_lyot_stop(
    wf::Proper.WaveFront,
    pupil_diameter_m::Real;
    pinhole_diameter_m::Real=pupil_diameter_m / 20,
    separation_m::Real=1.6 * pupil_diameter_m,
    position_angle_rad::Real=-pi / 2,
    obscuration_diameter_m::Real=0,
    center::Tuple{<:Real,<:Real}=(0, 0),
    outer_margin::Real=0.95,
    obscuration_margin::Real=1.15,
    spider_margin::Real=2,
)
    pupil_diameter_m > 0 || throw(ArgumentError("pupil_diameter_m must be positive"))
    pinhole_diameter_m >= 0 || throw(ArgumentError("pinhole_diameter_m must be nonnegative"))
    xc, yc = center
    main = prop_ellipse(
        wf,
        outer_margin * pupil_diameter_m / 2,
        outer_margin * pupil_diameter_m / 2,
        xc,
        yc,
    )
    if obscuration_diameter_m > 0
        main .*= prop_ellipse(
            wf,
            obscuration_margin * obscuration_diameter_m / 2,
            obscuration_margin * obscuration_diameter_m / 2,
            xc,
            yc;
            DARK=true,
        )
    end
    for polygon in subaru_spider_polygons(pupil_diameter_m / 2; margin=spider_margin)
        main .*= prop_irregular_polygon(
            wf,
            polygon.x .+ xc,
            polygon.y .+ yc;
            DARK=true,
        )
    end

    pinhole = zeros(eltype(main), size(main))
    if pinhole_diameter_m > 0
        pinhole_x = separation_m * cos(position_angle_rad) + xc
        pinhole_y = separation_m * sin(position_angle_rad) + yc
        pinhole = prop_ellipse(
            wf,
            pinhole_diameter_m / 2,
            pinhole_diameter_m / 2,
            pinhole_x,
            pinhole_y,
        )
    end

    intensity = abs2.(centered_field(wf))
    incident = sum(intensity)
    incident > 0 || throw(DomainError(incident, "Lyot plane has no incident flux"))
    amplitude_transmission = clamp.(main .+ pinhole, 0, 1)
    amplitude_reflection = sqrt.(max.(1 .- abs2.(amplitude_transmission), 0))
    return LyotStop(
        amplitude_transmission,
        amplitude_reflection,
        sum(abs2.(main) .* intensity) / incident,
        sum(abs2.(pinhole) .* intensity) / incident,
    )
end
