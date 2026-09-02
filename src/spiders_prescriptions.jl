const PROVISIONAL_SPIDERS_LLOWFS_SCHEMA =
    "org.subaru.spiders.llowfs-relative-intensity/provisional-1"
const PROVISIONAL_SPIDERS_SCC_SCHEMA =
    "org.subaru.spiders.scc-relative-intensity/provisional-1"

"""
Prepared numerical parameters for the provisional LLOWFS and SCC prescription.

`camera_magnification` is the detector resampling factor when one external
pupil sample maps to one propagation-grid sample. Execution compensates this
factor for any additional propagation-grid sampling so the caller-visible
detector coordinates do not change with `resolution`.
"""
struct SpidersSurrogateParameters{T<:Union{Float32,Float64}}
    beam_diameter_fraction::T
    camera_magnification::T
    detector_defocus_m::T
    central_tilt_cycles::T
    central_gaussian_phase_rad::T
    central_gaussian_sigma_ratio::T
end

struct ProvisionalSpidersLLOWFSPrescription{
    T<:Union{Float32,Float64},
} <: Function
    profile::ProvisionalSpidersProfile{T}
    parameters::SpidersSurrogateParameters{T}
end

struct ProvisionalSpidersSCCPrescription{
    T<:Union{Float32,Float64},
} <: Function
    profile::ProvisionalSpidersProfile{T}
    parameters::SpidersSurrogateParameters{T}
end

function _spiders_surrogate_parameters(
    ::Type{T};
    beam_diameter_fraction::Real,
    camera_magnification::Real,
    detector_defocus_m::Real,
    central_tilt_cycles::Real,
    central_gaussian_phase_rad::Real,
    central_gaussian_sigma_ratio::Real,
) where {T<:Union{Float32,Float64}}
    beam_fraction = T(beam_diameter_fraction)
    isfinite(beam_fraction) && zero(T) < beam_fraction < one(T) ||
        throw(ArgumentError(
            "beam_diameter_fraction must be finite and between zero and one",
        ))
    magnification = T(camera_magnification)
    isfinite(magnification) && magnification > zero(T) ||
        throw(ArgumentError(
            "camera_magnification must be finite and positive",
        ))
    defocus = T(detector_defocus_m)
    isfinite(defocus) || throw(ArgumentError(
        "detector_defocus_m must be finite",
    ))
    tilt = T(central_tilt_cycles)
    isfinite(tilt) || throw(ArgumentError(
        "central_tilt_cycles must be finite",
    ))
    gaussian = T(central_gaussian_phase_rad)
    isfinite(gaussian) || throw(ArgumentError(
        "central_gaussian_phase_rad must be finite",
    ))
    sigma = T(central_gaussian_sigma_ratio)
    isfinite(sigma) && sigma > zero(T) || throw(ArgumentError(
        "central_gaussian_sigma_ratio must be finite and positive",
    ))
    return SpidersSurrogateParameters{T}(
        beam_fraction,
        magnification,
        defocus,
        tilt,
        gaussian,
        sigma,
    )
end

function _validate_spiders_surrogate_grid(
    profile::ProvisionalSpidersProfile{T},
    parameters::SpidersSurrogateParameters{T},
    resolution::Integer,
    pupil_resolution::Integer,
    output_shape::Tuple{Int,Int},
) where {T}
    n = Int(resolution)
    n > 0 && iseven(n) || throw(ArgumentError(
        "the provisional SPIDERS propagation resolution must be positive and even",
    ))
    maximum(output_shape) <= n || throw(ArgumentError(
        "the provisional detector product $output_shape must fit inside the " *
        "$n-by-$n propagation grid",
    ))
    pupil_n = Int(pupil_resolution)
    pupil_n > 0 || throw(ArgumentError(
        "the provisional SPIDERS pupil resolution must be positive",
    ))

    optics = profile.optics
    pinhole_radius_norm = optics.reference_pinhole_diameter_m /
        optics.lyot_pupil_diameter_m
    pinhole_center_norm = T(2) *
        optics.reference_pinhole_separation_pupils
    grid_radius_norm = inv(parameters.beam_diameter_fraction)
    pinhole_center_norm + pinhole_radius_norm < grid_radius_norm ||
        throw(ArgumentError(
            "beam_diameter_fraction does not leave room for the SCC reference hole",
        ))
    pinhole_diameter_pixels = T(n) * parameters.beam_diameter_fraction *
        optics.reference_pinhole_diameter_m /
        optics.lyot_pupil_diameter_m
    pinhole_diameter_pixels >= T(2) || throw(ArgumentError(
        "the propagation grid must sample the SCC reference hole with at least two pixels",
    ))
    return n, pupil_n
end

function _spiders_pupil_resolution(
    resolution::Integer,
    pupil_resolution::Union{Nothing,Integer},
    beam_diameter_fraction::Real,
)
    if isnothing(pupil_resolution)
        return round(Int, resolution * beam_diameter_fraction)
    end
    return Int(pupil_resolution)
end

"""
    provisional_spiders_llowfs_configuration(profile; kwargs...)

Create an executable, explicitly unqualified LLOWFS Proper configuration. It
models the light rejected by a surrogate Lyot stop and then applies the
profile's provisional LLOWFS selector state. The as-built masks and detector
mapping can replace the prepared assets without changing the graph topology.
"""
function provisional_spiders_llowfs_configuration(
    profile::ProvisionalSpidersProfile{T};
    resolution::Integer=512,
    pupil_resolution::Union{Nothing,Integer}=nothing,
    beam_diameter_fraction::Real=0.25,
    camera_magnification::Real=0.10,
    detector_defocus_m::Real=0,
    central_tilt_cycles::Real=1,
    central_gaussian_phase_rad::Real=pi,
    central_gaussian_sigma_ratio::Real=0.35,
    rng_seed::Integer=0,
    output_schema::AbstractString=PROVISIONAL_SPIDERS_LLOWFS_SCHEMA,
) where {T}
    parameters = _spiders_surrogate_parameters(
        T;
        beam_diameter_fraction,
        camera_magnification,
        detector_defocus_m,
        central_tilt_cycles,
        central_gaussian_phase_rad,
        central_gaussian_sigma_ratio,
    )
    pupil_n = _spiders_pupil_resolution(
        resolution,
        pupil_resolution,
        beam_diameter_fraction,
    )
    n, pupil_n = _validate_spiders_surrogate_grid(
        profile,
        parameters,
        resolution,
        pupil_n,
        profile.llowfs_output_shape,
    )
    prescription = ProvisionalSpidersLLOWFSPrescription(
        profile,
        parameters,
    )
    return ProperPropagationConfiguration(
        prescription;
        resolution=n,
        pupil_resolution=pupil_n,
        output_rows=profile.llowfs_output_shape[1],
        output_columns=profile.llowfs_output_shape[2],
        diameter_m=profile.optics.apodizer_pupil_diameter_m,
        wavelength_um=profile.wavelength_um,
        rng_seed,
        output_schema,
    )
end

"""
    provisional_spiders_scc_configuration(profile; kwargs...)

Create an executable, explicitly unqualified SCC Proper configuration. The
captured profile selects main- and reference-beam transmission at preparation;
use separate prepared configurations for fringed and unfringed acquisition
states until a qualified chopper-state protocol is available.
"""
function provisional_spiders_scc_configuration(
    profile::ProvisionalSpidersProfile{T};
    resolution::Integer=512,
    pupil_resolution::Union{Nothing,Integer}=nothing,
    beam_diameter_fraction::Real=0.25,
    camera_magnification::Real=1,
    detector_defocus_m::Real=0,
    central_tilt_cycles::Real=1,
    central_gaussian_phase_rad::Real=pi,
    central_gaussian_sigma_ratio::Real=0.35,
    rng_seed::Integer=0,
    output_schema::AbstractString=PROVISIONAL_SPIDERS_SCC_SCHEMA,
) where {T}
    parameters = _spiders_surrogate_parameters(
        T;
        beam_diameter_fraction,
        camera_magnification,
        detector_defocus_m,
        central_tilt_cycles,
        central_gaussian_phase_rad,
        central_gaussian_sigma_ratio,
    )
    pupil_n = _spiders_pupil_resolution(
        resolution,
        pupil_resolution,
        beam_diameter_fraction,
    )
    n, pupil_n = _validate_spiders_surrogate_grid(
        profile,
        parameters,
        resolution,
        pupil_n,
        profile.scc_output_shape,
    )
    prescription = ProvisionalSpidersSCCPrescription(profile, parameters)
    return ProperPropagationConfiguration(
        prescription;
        resolution=n,
        pupil_resolution=pupil_n,
        output_rows=profile.scc_output_shape[1],
        output_columns=profile.scc_output_shape[2],
        diameter_m=profile.optics.apodizer_pupil_diameter_m,
        wavelength_um=profile.wavelength_um,
        rng_seed,
        output_schema,
    )
end

function _spiders_prescription_claims(
    configuration::ProperPropagationConfiguration,
    prescription,
    detector_choice::String,
)
    parameters = prescription.parameters
    profile = prescription.profile
    pupil_magnification = parameters.beam_diameter_fraction *
        configuration.resolution / configuration.pupil_resolution
    return (
        SpidersProfileClaim(
            :pupil_embedding,
            SpidersPlaceholder,
            "$(configuration.pupil_resolution)-pixel pupil embedded in a " *
            "$(configuration.resolution)-pixel grid at " *
            "$(parameters.beam_diameter_fraction) beam fraction " *
            "(resampling magnification $pupil_magnification)",
            "FFT padding is required to contain the off-axis SCC reference " *
            "hole; no qualified SPIDERS propagation-grid sampling was found",
            "Confirm the required pupil sampling, FFT padding, centering, and " *
            "anti-aliasing policy.",
        ),
        SpidersProfileClaim(
            :tilt_gaussian_vortex_surrogate,
            SpidersPlaceholder,
            "charge $(profile.focal_plane_mask.vortex_charge), " *
            "$(parameters.central_tilt_cycles) central tilt cycles, " *
            "$(parameters.central_gaussian_phase_rad) rad Gaussian phase, " *
            "sigma ratio $(parameters.central_gaussian_sigma_ratio)",
            "Only the mask family, charge, and 5.7 lambda/D diameter are " *
            "available; the central complex transmission is a numerical guess",
            "Replace the surrogate with the as-built complex mask map and " *
            "measured registration.",
        ),
        SpidersProfileClaim(
            :proper_relay,
            SpidersPlaceholder,
            "f/$(profile.optics.focal_plane_f_number) first focus; " *
            "$(profile.optics.fpm_to_l2_m) m FPM-to-L2; " *
            "$(profile.optics.l2_focal_length_m) m L2 and camera relay",
            "The optical design provides these headline dimensions, but not " *
            "a complete signed, folded, as-built surface prescription",
            "Provide the surface order, signed separations, conjugates, folds, " *
            "and measured aberration maps.",
        ),
        SpidersProfileClaim(
            :provisional_detector_model,
            SpidersPlaceholder,
            detector_choice,
            "Local calibration products establish legacy array shapes only; " *
            "the surrogate emits relative monochromatic intensity",
            "Provide the crop, plate scale, registration, throughput, exposure, " *
            "gain, noise, nonlinearity, bad-pixel, and saturation models.",
        ),
    )
end

"""Return the numerical assumptions made by one provisional prescription."""
function spiders_prescription_claims(
    configuration::ProperPropagationConfiguration{
        T,
        <:ProvisionalSpidersLLOWFSPrescription,
    },
) where {T}
    prescription = configuration.prescription
    detector_choice = "LLOWFS $(configuration.output_rows)-by-" *
        "$(configuration.output_columns), $(prescription.profile.llowfs_plane), " *
        "$(prescription.profile.llowfs_magnification), magnification " *
        "$(prescription.parameters.camera_magnification), defocus " *
        "$(prescription.parameters.detector_defocus_m) m"
    return _spiders_prescription_claims(
        configuration,
        prescription,
        detector_choice,
    )
end

function spiders_prescription_claims(
    configuration::ProperPropagationConfiguration{
        T,
        <:ProvisionalSpidersSCCPrescription,
    },
) where {T}
    prescription = configuration.prescription
    detector_choice = "SCC $(configuration.output_rows)-by-" *
        "$(configuration.output_columns), $(prescription.profile.scc_plane), " *
        "magnification $(prescription.parameters.camera_magnification), " *
        "defocus $(prescription.parameters.detector_defocus_m) m, rotation " *
        "$(prescription.profile.optics.camera_rotation_deg) deg"
    return _spiders_prescription_claims(
        configuration,
        prescription,
        detector_choice,
    )
end

"""Return the complete profile and numerical qualification checklist."""
function spiders_configuration_claims(
    configuration::ProperPropagationConfiguration{
        T,
        <:Union{
            ProvisionalSpidersLLOWFSPrescription,
            ProvisionalSpidersSCCPrescription,
        },
    },
) where {T}
    return (
        spiders_profile_claims(configuration.prescription.profile)...,
        spiders_prescription_claims(configuration)...,
    )
end

@inline function _spiders_centered_coordinate(
    index::Int,
    resolution::Int,
    beam_radius_pixels::T,
) where {T}
    return (T(index) - T(resolution ÷ 2 + 1)) / beam_radius_pixels
end

function _spiders_surrogate_masks(
    profile::ProvisionalSpidersProfile{T},
    parameters::SpidersSurrogateParameters{T},
    resolution::Int,
) where {T}
    optics = profile.optics
    beam_radius_pixels = T(resolution) *
        parameters.beam_diameter_fraction / T(2)
    focal_radius_pixels = profile.focal_plane_mask.diameter_lambda_over_d /
        (T(2) * parameters.beam_diameter_fraction)
    reference_radius_norm = optics.reference_pinhole_diameter_m /
        optics.lyot_pupil_diameter_m
    reference_offset_norm = T(2) *
        optics.reference_pinhole_separation_pupils
    reference_angle = deg2rad(optics.reference_pinhole_position_angle_deg)
    reference_x = reference_offset_norm * cos(reference_angle)
    reference_y = reference_offset_norm * sin(reference_angle)

    apodizer = Matrix{T}(undef, resolution, resolution)
    focal_plane_mask = Matrix{Complex{T}}(
        undef,
        resolution,
        resolution,
    )
    main_lyot = Matrix{T}(undef, resolution, resolution)
    reference_hole = Matrix{T}(undef, resolution, resolution)

    @inbounds for column in 1:resolution
        x_norm = _spiders_centered_coordinate(
            column,
            resolution,
            beam_radius_pixels,
        )
        x_focal = T(column) - T(resolution ÷ 2 + 1)
        for row in 1:resolution
            y_norm = _spiders_centered_coordinate(
                row,
                resolution,
                beam_radius_pixels,
            )
            y_focal = T(row) - T(resolution ÷ 2 + 1)
            pupil_radius = hypot(x_norm, y_norm)
            apodizer[row, column] = T(
                optics.apodizer_inner_diameter_ratio <= pupil_radius <=
                optics.apodizer_outer_diameter_ratio,
            )

            focal_radius = hypot(x_focal, y_focal)
            focal_angle = atan(y_focal, x_focal)
            phase = T(profile.focal_plane_mask.vortex_charge) * focal_angle
            if focal_radius <= focal_radius_pixels
                normalized_x = x_focal / (T(2) * focal_radius_pixels)
                normalized_radius = focal_radius / focal_radius_pixels
                phase += T(2pi) * parameters.central_tilt_cycles *
                    normalized_x
                phase += parameters.central_gaussian_phase_rad * exp(
                    -T(0.5) * (
                        normalized_radius /
                        parameters.central_gaussian_sigma_ratio
                    )^2,
                )
            end
            focal_plane_mask[row, column] = cis(phase)

            lyot_spider = (
                abs(x_norm) <= optics.lyot_spider_width_ratio ||
                abs(y_norm) <= optics.lyot_spider_width_ratio
            )
            main_lyot[row, column] = T(
                optics.lyot_inner_diameter_ratio <= pupil_radius <=
                optics.lyot_outer_diameter_ratio && !lyot_spider,
            )
            reference_hole[row, column] = T(
                hypot(
                    x_norm - reference_x,
                    y_norm - reference_y,
                ) <= reference_radius_norm,
            )
        end
    end
    return (; apodizer, focal_plane_mask, main_lyot, reference_hole)
end

function _spiders_lyot_mask(
    prescription::ProvisionalSpidersLLOWFSPrescription{T},
    main_lyot::Matrix{T},
    reference_hole::Matrix{T},
) where {T}
    lyot = similar(main_lyot)
    @. lyot = max(zero(T), one(T) - main_lyot - reference_hole)
    return lyot
end

function _spiders_lyot_mask(
    prescription::ProvisionalSpidersSCCPrescription{T},
    main_lyot::Matrix{T},
    reference_hole::Matrix{T},
) where {T}
    main_scale = prescription.profile.scc_main_beam === :unblocked ?
        one(T) : zero(T)
    reference_scale = prescription.profile.scc_reference_beam === :open ?
        one(T) : zero(T)
    lyot = similar(main_lyot)
    @. lyot = main_scale * main_lyot + reference_scale * reference_hole
    return lyot
end

function _spiders_target_array(target, values::Matrix{T}) where {T}
    prepared = _Backends.allocate_device_array(
        target,
        T,
        size(values)...,
    )
    copyto!(prepared, values)
    return prepared
end

function _prepare_spiders_common_assets(
    profile::ProvisionalSpidersProfile{T},
    parameters::SpidersSurrogateParameters{T},
    target,
    resolution::Int,
) where {T}
    masks = _spiders_surrogate_masks(profile, parameters, resolution)
    focal_plane_mask_internal = Proper.prop_shift_center(
        masks.focal_plane_mask;
        inverse=true,
    )
    square_intensity = _Backends.allocate_device_array(
        target,
        T,
        resolution,
        resolution,
    )
    rotated_intensity = similar(square_intensity)
    centered_field = _Backends.allocate_device_array(
        target,
        Complex{T},
        resolution,
        resolution,
    )
    padded_pupil_opd = similar(square_intensity)
    padded_pupil_amplitude = similar(square_intensity)
    fill!(square_intensity, zero(T))
    fill!(rotated_intensity, zero(T))
    fill!(centered_field, zero(Complex{T}))
    fill!(padded_pupil_opd, zero(T))
    fill!(padded_pupil_amplitude, zero(T))
    assets = (
        apodizer_mask=_spiders_target_array(target, masks.apodizer),
        focal_plane_mask_internal=_spiders_target_array(
            target,
            focal_plane_mask_internal,
        ),
        padded_pupil_opd,
        padded_pupil_amplitude,
        square_intensity,
        rotated_intensity,
        centered_field,
        pupil_resample_context=Proper.RunContext(typeof(square_intensity)),
        detector_resample_context=Proper.RunContext(typeof(square_intensity)),
    )
    return masks, assets
end

@inline _spiders_output_shape(
    prescription::ProvisionalSpidersLLOWFSPrescription,
) = prescription.profile.llowfs_output_shape

@inline _spiders_output_shape(
    prescription::ProvisionalSpidersSCCPrescription,
) = prescription.profile.scc_output_shape

function prepare_proper_assets(
    prescription::Union{
        ProvisionalSpidersLLOWFSPrescription{T},
        ProvisionalSpidersSCCPrescription{T},
    },
    target,
    ::Type{T},
    resolution::Int,
    output_shape::Tuple{Int,Int},
) where {T}
    expected_shape = _spiders_output_shape(prescription)
    output_shape == expected_shape || throw(DimensionMismatch(
        "SPIDERS output shape $output_shape does not match $expected_shape",
    ))
    masks, assets = _prepare_spiders_common_assets(
        prescription.profile,
        prescription.parameters,
        target,
        resolution,
    )
    lyot = _spiders_lyot_mask(
        prescription,
        masks.main_lyot,
        masks.reference_hole,
    )
    return merge(assets, (
        lyot_mask=_spiders_target_array(target, lyot),
    ))
end

@inline _spiders_rotate_camera(
    ::ProvisionalSpidersLLOWFSPrescription,
    square_intensity,
    rotated_intensity,
    run_context,
) = square_intensity

@inline function _spiders_rotate_scc_camera(
    profile::ProvisionalSpidersProfile,
    square_intensity,
    rotated_intensity,
    run_context,
)
    Proper.prop_rotate!(
        rotated_intensity,
        square_intensity,
        profile.optics.camera_rotation_deg,
        run_context,
    )
    return rotated_intensity
end

@inline _spiders_rotate_camera(
    prescription::ProvisionalSpidersSCCPrescription,
    square_intensity,
    rotated_intensity,
    run_context,
) = _spiders_rotate_scc_camera(
    prescription.profile,
    square_intensity,
    rotated_intensity,
    run_context,
)

@inline function _spiders_camera_plane!(
    prescription::ProvisionalSpidersLLOWFSPrescription,
    wavefront,
    run_context,
)
    prescription.profile.llowfs_plane === :pupil_plane && return wavefront
    focal_length = prescription.profile.optics.l2_focal_length_m
    _prop_lens!(wavefront, focal_length, run_context, "LLOWFS camera relay")
    _prop_propagate!(
        wavefront,
        focal_length + prescription.parameters.detector_defocus_m,
        run_context,
        "LLOWFS detector",
    )
    return wavefront
end

@inline function _spiders_scc_camera_plane!(
    profile::ProvisionalSpidersProfile,
    parameters::SpidersSurrogateParameters,
    wavefront,
    run_context,
)
    profile.scc_plane === :pupil_plane && return wavefront
    focal_length = profile.optics.l2_focal_length_m
    _prop_lens!(wavefront, focal_length, run_context, "SCC camera relay")
    _prop_propagate!(
        wavefront,
        focal_length + parameters.detector_defocus_m,
        run_context,
        "SCC detector",
    )
    return wavefront
end

@inline _spiders_camera_plane!(
    prescription::ProvisionalSpidersSCCPrescription,
    wavefront,
    run_context,
) = _spiders_scc_camera_plane!(
    prescription.profile,
    prescription.parameters,
    wavefront,
    run_context,
)

function _spiders_propagate_to_lyot!(
    prescription,
    wavelength_m,
    resolution;
    pupil_opd,
    pupil_amplitude,
    diameter_m,
    field,
    wavefront,
    run_context,
    apodizer_mask,
    focal_plane_mask_internal,
    padded_pupil_opd,
    padded_pupil_amplitude,
    pupil_resample_context,
)
    size(field) == (resolution, resolution) || throw(DimensionMismatch(
        "SPIDERS field binding changed after preparation",
    ))
    parameters = prescription.parameters
    optics = prescription.profile.optics
    pupil_magnification = parameters.beam_diameter_fraction *
        typeof(parameters.beam_diameter_fraction)(resolution) /
        typeof(parameters.beam_diameter_fraction)(size(pupil_opd, 1))
    Proper.prop_magnify!(
        padded_pupil_opd,
        pupil_opd,
        pupil_magnification,
        pupil_resample_context,
    )
    Proper.prop_magnify!(
        padded_pupil_amplitude,
        pupil_amplitude,
        pupil_magnification,
        pupil_resample_context,
    )
    Proper.prop_begin!(
        wavefront,
        diameter_m,
        wavelength_m;
        beam_diam_fraction=parameters.beam_diameter_fraction,
    )
    Proper.prop_multiply(wavefront, padded_pupil_amplitude)
    Proper.prop_add_phase(wavefront, padded_pupil_opd)
    Proper.prop_multiply(wavefront, apodizer_mask)

    focal_length = optics.focal_plane_f_number * diameter_m
    _prop_lens!(
        wavefront,
        focal_length,
        run_context,
        "SPIDERS focal-plane relay",
    )
    _prop_propagate!(
        wavefront,
        focal_length,
        run_context,
        "SPIDERS focal-plane mask",
    )
    wavefront.field .*= focal_plane_mask_internal

    _prop_propagate!(
        wavefront,
        optics.fpm_to_l2_m,
        run_context,
        "SPIDERS L2",
    )
    _prop_lens!(
        wavefront,
        optics.l2_focal_length_m,
        run_context,
        "SPIDERS L2",
    )
    _prop_propagate!(
        wavefront,
        optics.l2_focal_length_m,
        run_context,
        "SPIDERS Lyot stop",
    )
    return pupil_magnification
end

function _spiders_propagate_from_lyot!(
    prescription,
    pupil_magnification,
    wavefront;
    output,
    run_context,
    lyot_mask,
    square_intensity,
    rotated_intensity,
    centered_field,
    detector_resample_context,
)
    parameters = prescription.parameters
    Proper.prop_multiply(wavefront, lyot_mask)
    _spiders_camera_plane!(prescription, wavefront, run_context)

    Proper.prop_shift_center!(centered_field, wavefront.field)
    @. square_intensity = abs2(centered_field)
    sampling_m = Proper.prop_get_sampling(wavefront)
    camera_intensity = _spiders_rotate_camera(
        prescription,
        square_intensity,
        rotated_intensity,
        run_context,
    )
    detector_magnification =
        parameters.camera_magnification / pupil_magnification
    Proper.prop_magnify!(
        output,
        camera_intensity,
        detector_magnification,
        detector_resample_context,
    )
    return output, sampling_m / detector_magnification
end

function _run_provisional_spiders_prescription(
    prescription,
    wavelength_m,
    resolution;
    pupil_opd,
    pupil_amplitude,
    diameter_m,
    field,
    wavefront,
    output,
    run_context,
    apodizer_mask,
    focal_plane_mask_internal,
    lyot_mask,
    padded_pupil_opd,
    padded_pupil_amplitude,
    square_intensity,
    rotated_intensity,
    centered_field,
    pupil_resample_context,
    detector_resample_context,
)
    pupil_magnification = _spiders_propagate_to_lyot!(
        prescription,
        wavelength_m,
        resolution;
        pupil_opd,
        pupil_amplitude,
        diameter_m,
        field,
        wavefront,
        run_context,
        apodizer_mask,
        focal_plane_mask_internal,
        padded_pupil_opd,
        padded_pupil_amplitude,
        pupil_resample_context,
    )
    return _spiders_propagate_from_lyot!(
        prescription,
        pupil_magnification,
        wavefront;
        output,
        run_context,
        lyot_mask,
        square_intensity,
        rotated_intensity,
        centered_field,
        detector_resample_context,
    )
end

function (prescription::ProvisionalSpidersLLOWFSPrescription)(
    wavelength_m,
    resolution;
    kwargs...,
)
    return _run_provisional_spiders_prescription(
        prescription,
        wavelength_m,
        resolution;
        kwargs...,
    )
end

function (prescription::ProvisionalSpidersSCCPrescription)(
    wavelength_m,
    resolution;
    kwargs...,
)
    return _run_provisional_spiders_prescription(
        prescription,
        wavelength_m,
        resolution;
        kwargs...,
    )
end
