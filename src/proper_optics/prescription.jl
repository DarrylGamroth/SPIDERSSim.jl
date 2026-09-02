"""
    SCCPropagationResult

Contain the final PROPER wavefront, Lyot-stop transmission metrics, and the
fractional throughput retained by each finite optic aperture.
"""
struct SCCPropagationResult{W,L}
    wavefront::W
    lyot_stop::L
    aperture_throughput::NamedTuple{(:oae1, :oae2, :lens0, :lens1),NTuple{4,Float64}}
end

"""Return the centered final complex electric field from a SPIDERS result."""
scc_field(result::SCCPropagationResult) = centered_field(result.wavefront)

"""Return the centered final intensity from a SPIDERS result."""
scc_intensity(result::SCCPropagationResult) = abs2.(scc_field(result))

# Proper.jl's registered 0.1.0 accessor returned internal FFT ordering while
# the current checkout returns a centered copy. Build the application-level
# contract from the stable public centering primitive so both revisions agree.
centered_field(wf::Proper.WaveFront) = prop_shift_center(wf.field)

function lens_aperture!(wf::Proper.WaveFront, focal_length_m::Real, diameter_m::Real, name::AbstractString)
    input_power = sum(abs2, wf.field)
    prop_lens(wf, focal_length_m, name)
    prop_circular_aperture(wf, diameter_m / 2)
    output_power = sum(abs2, wf.field)
    return Float64(output_power / input_power)
end

function apply_entrance_pupil!(wf::Proper.WaveFront, configuration::SCCPropagationConfiguration)
    diameter = configuration.telescope_diameter_m
    prop_circular_aperture(wf, diameter / 2)
    prop_define_entrance(wf)
    if configuration.pupil_mode === :fits
        prop_errormap(
            wf,
            configuration.pupil_path;
            AMPLITUDE=true,
            SAMPLING=diameter / 1200,
        )
    else
        subaru_pupil!(
            wf;
            obscuration_ratio=configuration.pupil_obscuration_ratio,
            include_spiders=configuration.include_spiders,
        )
    end
    return wf
end

function apply_ao_residual!(
    wf::Proper.WaveFront,
    configuration::SCCPropagationConfiguration,
    ao_residual_nm::Union{Nothing,AbstractMatrix},
    ao_sampling_m::Union{Nothing,Real},
)
    if ao_residual_nm !== nothing
        sampling = something(ao_sampling_m, 0.045 / 7.73 * configuration.telescope_diameter_m)
        center_x = size(ao_residual_nm, 2) / 2
        center_y = size(ao_residual_nm, 1) / 2
        residual = prop_resamplemap(wf, ao_residual_nm, sampling, center_x, center_y)
        prop_add_phase(wf, residual .* 1e-9)
    elseif configuration.turbulence
        prop_errormap(
            wf,
            configuration.ao_residual_path;
            WAVEFRONT=true,
            SAMPLING=0.045 / 7.73 * configuration.telescope_diameter_m,
            NM=true,
        )
    end
    return wf
end

function apply_apodizer!(wf::Proper.WaveFront, configuration::SCCPropagationConfiguration)
    configuration.coronagraph || return wf
    configuration.apodizer_mode === :none && return wf
    if configuration.apodizer_mode === :fits
        prop_errormap(
            wf,
            configuration.apodizer_path;
            AMPLITUDE=true,
            SAMPLING=2 * prop_get_beamradius(wf) / 1200,
        )
    else
        apodizer = radial_apodizer(
            configuration.apodizer_polynomial,
            prop_get_gridsize(wf),
            prop_get_beamradius(wf),
            configuration.pupil_obscuration_ratio,
            prop_get_sampling(wf),
        )
        prop_multiply(wf, apodizer)
    end
    return wf
end

function apply_fpm!(wf::Proper.WaveFront, configuration::SCCPropagationConfiguration)
    configuration.coronagraph || return wf
    wavelength = prop_get_wavelength(wf)
    if configuration.fpm_mode === :analytic
        radius = configuration.fpm_band === :J ? 459.62e-6 / 2 : 584.63e-6 / 2
        if configuration.fpm_type === :TG
            phase = tilt_gaussian_fpm(
                wf;
                radius_m=radius,
                tilt_rad=1.6 * atan(1 / 64),
                gaussian_phase=7.8,
                charge=0,
                piston_rad=1.9999pi,
                levels=16,
                center_pupils=configuration.center_pupils,
            )
        else
            sampling = prop_get_sampling(wf)
            angular_sampling = prop_get_sampling_radians(wf)
            resolution_element_m = wavelength / wf.beam_diameter_m / angular_sampling * sampling
            phase = tilt_gaussian_fpm(
                wf;
                radius_m=resolution_element_m,
                tilt_rad=1.6 * atan(1 / 64),
                gaussian_phase=10,
                charge=4,
                piston_rad=0,
                levels=16,
                center_pupils=configuration.center_pupils,
            )
        end
        prop_add_phase(wf, phase .* (wavelength / 2pi))
    elseif configuration.fpm_map_type === :mirror_surface
        prop_errormap(
            wf,
            something(configuration.fpm_map_path);
            MIRROR_SURFACE=true,
            SAMPLING=configuration.fpm_map_sampling_m,
            ROTATEMAP=180,
        )
    else
        multiplier = configuration.fpm_map_type === :phase ? wavelength / 2pi : 1.0
        prop_errormap(
            wf,
            something(configuration.fpm_map_path);
            WAVEFRONT=true,
            SAMPLING=configuration.fpm_map_sampling_m,
            MULTIPLY=multiplier,
        )
    end

    if configuration.fpm_mode === :map && configuration.center_pupils
        n = prop_get_gridsize(wf)
        axis = ((0:(n - 1)) .- n / 2) .* prop_get_sampling(wf)
        yy = reshape(axis, :, 1)
        tilt = repeat((-2pi * 1.6 * atan(1 / 64) / wavelength) .* yy, 1, n)
        prop_add_phase(wf, -tilt .* (wavelength / (4pi)))
    end
    return wf
end

"""
    propagate_scc(wavelength_m, gridsize; configuration=SCCPropagationConfiguration(),
                   ao_residual_nm=nothing, ao_sampling_m=nothing, rng=default_rng())

Propagate a monochromatic complex field through the SPIDERS SCC prescription
and return a [`SCCPropagationResult`](@ref).

`ao_residual_nm` may supply one in-memory AO residual frame. Its physical pixel
spacing defaults to the PASSATA/GPI scale used by the MATLAB movie model and
can be overridden with `ao_sampling_m`.
"""
function propagate_scc(
    wavelength_m::Real,
    gridsize::Integer;
    configuration::SCCPropagationConfiguration=SCCPropagationConfiguration(),
    ao_residual_nm::Union{Nothing,AbstractMatrix}=nothing,
    ao_sampling_m::Union{Nothing,Real}=nothing,
    rng::AbstractRNG=default_rng(),
)
    validate(configuration)
    wavelength_m > 0 || throw(ArgumentError("wavelength_m must be positive"))
    gridsize > 0 || throw(ArgumentError("gridsize must be positive"))

    diameter = configuration.telescope_diameter_m
    input_focal_length = diameter * configuration.input_fratio
    wf = prop_begin(
        diameter,
        wavelength_m,
        gridsize;
        beam_diam_fraction=configuration.beam_diameter_fraction,
    )
    apply_entrance_pupil!(wf, configuration)
    apply_ao_residual!(wf, configuration, ao_residual_nm, ao_sampling_m)
    if configuration.scintillation
        prop_psd_errormap(
            wf,
            0.01,
            22 / diameter,
            11 / 3;
            AMPLITUDE=1,
            TPF=true,
            RNG=rng,
        )
    end

    prop_propagate(wf, input_focal_length, "dummy lens")
    prop_lens(wf, input_focal_length, "F/13.901 convergent beam")
    prop_propagate(wf, input_focal_length, "Bay #4 focus")

    prop_propagate(wf, 772.35e-3, "OAE1")
    oae1_throughput = lens_aperture!(wf, 459.7e-3, 3 * 25.4e-3, "OAE1")
    if configuration.oae1_error_path !== nothing
        prop_errormap(
            wf,
            configuration.oae1_error_path;
            MULTIPLY=1,
            MIRROR_SURFACE=true,
            SAMPLING=65.88e-6,
            MICRONS=true,
        )
    end

    prop_propagate(wf, 471.4e-3, "DM")
    if configuration.oae1_error_path !== nothing && configuration.dm_correction
        wavefront_error = angle.(centered_field(wf)) .* (prop_get_wavelength(wf) / 2pi)
        actuator_spacing = 2 * prop_get_beamradius(wf) / configuration.dm_actuators_across_pupil
        dm_map = prop_magnify(
            wavefront_error,
            prop_get_sampling(wf) / actuator_spacing,
            configuration.dm_actuators,
        )
        actuator_center = configuration.dm_actuators / 2
        prop_dm(
            wf,
            -dm_map ./ 2,
            actuator_center,
            actuator_center,
            actuator_spacing;
            FIT=true,
        )
    end

    distance_to_focus = prop_get_distancetofocus(wf)
    prop_propagate(wf, distance_to_focus + 367.2e-3, "OAE2")
    oae2_throughput = lens_aperture!(wf, 277.8e-3, 38.1e-3, "OAE2")
    if configuration.oae2_error_path !== nothing
        prop_errormap(
            wf,
            configuration.oae2_error_path;
            MULTIPLY=1,
            MIRROR_SURFACE=true,
            SAMPLING=126e-6,
            MICRONS=true,
        )
    end

    prop_propagate(wf, 377.2e-3, "Apodizer")
    apply_apodizer!(wf, configuration)
    prop_propagate(wf, prop_get_distancetofocus(wf), "FPM")
    apply_fpm!(wf, configuration)

    prop_propagate(wf, 358.69e-3 + 5.166e-3, "Lens 0")
    lens0_throughput = lens_aperture!(wf, 285.99e-3, 40e-3, "Lens 0")
    prop_propagate(wf, 382.67e-3 - 0.022e-3, "Lyot")

    pupil_diameter = configuration.lyot_pupil_diameter_m
    center_y = configuration.center_pupils ? configuration.lyot_reference_separation * pupil_diameter / 2 : 0.0
    pinhole_diameter = configuration.reference_pinhole ? configuration.lyot_pinhole_ratio * pupil_diameter : 0.0
    lyot = scc_lyot_stop(
        wf,
        pupil_diameter;
        pinhole_diameter_m=pinhole_diameter,
        separation_m=configuration.lyot_reference_separation * pupil_diameter,
        position_angle_rad=configuration.lyot_reference_angle_rad,
        obscuration_diameter_m=configuration.pupil_obscuration_ratio * pupil_diameter,
        center=(0.0, center_y),
        outer_margin=configuration.lyot_outer_margin,
        obscuration_margin=configuration.lyot_obscuration_margin,
        spider_margin=configuration.lyot_spider_margin,
    )
    configuration.coronagraph && prop_multiply(wf, lyot.amplitude_transmission)

    prop_propagate(wf, 55e-3, "Chopper")
    prop_propagate(wf, 133.16e-3 - 0.022e-3, "Lens 1")
    lens1_throughput = lens_aperture!(wf, 285.99e-3, 40e-3, "Lens 1")
    prop_propagate(wf, prop_get_distancetofocus(wf), "SCC Focus")

    throughput = (
        oae1=oae1_throughput,
        oae2=oae2_throughput,
        lens0=lens0_throughput,
        lens1=lens1_throughput,
    )
    return SCCPropagationResult(wf, lyot, throughput)
end

"""Compatibility alias for [`propagate_scc`](@ref)."""
