"""
    LLOWFSConfig(; kwargs...)

Configure the reflected-Lyot Low-Order Wavefront Sensor relay. The defaults
are the Zemax configuration-5 prescription, converted from millimetres to
metres. Propagation through glass uses reduced distances (`thickness/index`),
so the three paraxial refracting surfaces of the cemented doublet can be
represented with PROPER's vacuum-wavelength propagation.

The -5 degree coordinate break after the Lyot mirror only redirects the local
optical axis and is therefore not applied to the centered scalar field. The
Zemax model contains no measured reflective-mask coating or figure map.
"""
Base.@kwdef struct LLOWFSConfig
    temperature_c::Float64 = 20.0
    lyot_to_lens_m::Float64 = 515e-3
    front_curvature_inv_m::Float64 = 7.6753039420361041
    cement_curvature_inv_m::Float64 = -22.355360815523563
    s_lah71_thickness_m::Float64 = 5.22416e-3
    s_nph5_thickness_m::Float64 = 1.8e-3
    lens_clear_radius_m::Float64 = 8e-3
    lens_to_filter_m::Float64 = 150e-3
    filter_thickness_m::Float64 = 2e-3
    filter_to_detector_m::Float64 = 71.8759395036e-3
    detector_width_m::Float64 = 4.8e-3
    detector_height_m::Float64 = 3.84e-3
    detector_pixel_pitch_m::Float64 = 15e-6
    detector_pixels_x::Int = 320
    detector_pixels_y::Int = 256
    measured_roi_pixels_x::Int = 34
    measured_roi_pixels_y::Int = 47
end

function validate(config::LLOWFSConfig)
    config.lyot_to_lens_m >= 0 || throw(ArgumentError(
        "lyot_to_lens_m must be nonnegative"))
    config.s_lah71_thickness_m > 0 || throw(ArgumentError(
        "s_lah71_thickness_m must be positive"))
    config.s_nph5_thickness_m > 0 || throw(ArgumentError(
        "s_nph5_thickness_m must be positive"))
    config.lens_clear_radius_m > 0 || throw(ArgumentError(
        "lens_clear_radius_m must be positive"))
    config.lens_to_filter_m >= 0 || throw(ArgumentError(
        "lens_to_filter_m must be nonnegative"))
    config.filter_thickness_m >= 0 || throw(ArgumentError(
        "filter_thickness_m must be nonnegative"))
    config.filter_to_detector_m >= 0 || throw(ArgumentError(
        "filter_to_detector_m must be nonnegative"))
    config.detector_width_m > 0 || throw(ArgumentError(
        "detector_width_m must be positive"))
    config.detector_height_m > 0 || throw(ArgumentError(
        "detector_height_m must be positive"))
    config.detector_pixel_pitch_m > 0 || throw(ArgumentError(
        "detector_pixel_pitch_m must be positive"))
    config.detector_pixels_x > 0 || throw(ArgumentError(
        "detector_pixels_x must be positive"))
    config.detector_pixels_y > 0 || throw(ArgumentError(
        "detector_pixels_y must be positive"))
    config.measured_roi_pixels_x > 0 || throw(ArgumentError(
        "measured_roi_pixels_x must be positive"))
    config.measured_roi_pixels_y > 0 || throw(ArgumentError(
        "measured_roi_pixels_y must be positive"))
    isapprox(
        config.detector_width_m,
        config.detector_pixels_x * config.detector_pixel_pitch_m;
        rtol=0,
        atol=eps(config.detector_width_m) * 8,
    ) || throw(ArgumentError(
        "detector width, pixel count, and pixel pitch are inconsistent"))
    isapprox(
        config.detector_height_m,
        config.detector_pixels_y * config.detector_pixel_pitch_m;
        rtol=0,
        atol=eps(config.detector_height_m) * 8,
    ) || throw(ArgumentError(
        "detector height, pixel count, and pixel pitch are inconsistent"))
    return config
end

# Zemax Sellmeier-1 coefficients from the glass catalogs embedded in the ZAR.
# Coefficient order is B1, C1, B2, C2, B3, C3; wavelength is in micrometres.
const _S_LAH71_DISPERSION = (
    1.982800310,
    1.189874560e-2,
    3.167584500e-1,
    5.271560010e-2,
    2.444726460,
    2.132206970e2,
)
const _S_NPH5_DISPERSION = (
    1.891089960,
    1.411644990e-2,
    3.952201260e-1,
    6.628344450e-2,
    2.204921270,
    1.486807000e2,
)
const _F_SILICA_DISPERSION = (
    6.961663000e-1,
    4.679148000e-3,
    4.079426000e-1,
    1.351206300e-2,
    8.974794000e-1,
    9.793400250e1,
)

# Zemax thermal data: D0, D1, D2, E0, E1, lambda_TK, reference temperature.
const _S_LAH71_THERMAL = (
    -3.710000000e-7,
    9.010000000e-9,
    -2.990000000e-11,
    9.140000000e-7,
    1.110000000e-9,
    2.570000000e-1,
    25.0,
)
const _S_NPH5_THERMAL = (
    -3.540000000e-6,
    1.020000000e-8,
    -1.620000000e-11,
    9.630000000e-7,
    1.370000000e-9,
    3.120000000e-1,
    25.0,
)
const _F_SILICA_THERMAL = (
    2.237000000e-5,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    20.0,
)

@inline function _zemax_sellmeier1_index(
    wavelength_m::Real,
    temperature_c::Real,
    dispersion::NTuple{6,<:Real},
    thermal::NTuple{7,<:Real},
)
    wavelength_um = wavelength_m * 1e6
    wavelength_um > 0 || throw(ArgumentError("wavelength_m must be positive"))
    wavelength2 = wavelength_um * wavelength_um
    b1, c1, b2, c2, b3, c3 = dispersion
    index = sqrt(
        1 +
        b1 * wavelength2 / (wavelength2 - c1) +
        b2 * wavelength2 / (wavelength2 - c2) +
        b3 * wavelength2 / (wavelength2 - c3),
    )
    d0, d1, d2, e0, e1, lambda_tk, reference_temperature_c = thermal
    delta_temperature = temperature_c - reference_temperature_c
    thermal_term =
        d0 * delta_temperature +
        d1 * delta_temperature^2 +
        d2 * delta_temperature^3 +
        (e0 * delta_temperature + e1 * delta_temperature^2) /
        (wavelength2 - lambda_tk^2)
    return index + (index^2 - 1) / (2index) * thermal_term
end

"""
    llowfs_refractive_indices(wavelength_m; temperature_c=20)

Return the Zemax-catalog refractive indices for the two LLOWFS doublet glasses
and fused-silica filter at `wavelength_m`.
"""
function llowfs_refractive_indices(
    wavelength_m::Real;
    temperature_c::Real=20.0,
)
    return (
        s_lah71=_zemax_sellmeier1_index(
            wavelength_m,
            temperature_c,
            _S_LAH71_DISPERSION,
            _S_LAH71_THERMAL,
        ),
        s_nph5=_zemax_sellmeier1_index(
            wavelength_m,
            temperature_c,
            _S_NPH5_DISPERSION,
            _S_NPH5_THERMAL,
        ),
        fused_silica=_zemax_sellmeier1_index(
            wavelength_m,
            temperature_c,
            _F_SILICA_DISPERSION,
            _F_SILICA_THERMAL,
        ),
    )
end

@inline _translation_matrix(distance) = (one(distance), distance, zero(distance), one(distance))
@inline _refraction_matrix(power) = (one(power), zero(power), -power, one(power))

@inline function _multiply_abcd(left, right)
    a, b, c, d = left
    e, f, g, h = right
    return (
        muladd(a, e, b * g),
        muladd(a, f, b * h),
        muladd(c, e, d * g),
        muladd(c, f, d * h),
    )
end

function _llowfs_lens_matrix(wavelength_m::Real, config::LLOWFSConfig)
    indices = llowfs_refractive_indices(
        wavelength_m;
        temperature_c=config.temperature_c,
    )
    front_power = (indices.s_lah71 - 1) * config.front_curvature_inv_m
    cement_power = (indices.s_nph5 - indices.s_lah71) *
                   config.cement_curvature_inv_m
    return _multiply_abcd(
        _translation_matrix(config.s_nph5_thickness_m / indices.s_nph5),
        _multiply_abcd(
            _refraction_matrix(cement_power),
            _multiply_abcd(
                _translation_matrix(config.s_lah71_thickness_m / indices.s_lah71),
                _refraction_matrix(front_power),
            ),
        ),
    ), indices
end

@inline function _cement_to_detector_reduced_m(config::LLOWFSConfig, indices)
    return config.s_nph5_thickness_m / indices.s_nph5 +
           config.lens_to_filter_m +
           config.filter_thickness_m / indices.fused_silica +
           config.filter_to_detector_m
end

"""
    llowfs_relay_summary(wavelength_m; config=LLOWFSConfig())

Return first-order checks derived from the Zemax LLOWFS prescription. Distances
are measured from the final doublet surface. A negative `detector_defocus_m`
means the detector is before the paraxial conjugate.
"""
function llowfs_relay_summary(
    wavelength_m::Real;
    config::LLOWFSConfig=LLOWFSConfig(),
)
    validate(config)
    lens, indices = _llowfs_lens_matrix(wavelength_m, config)
    a, _, c, _ = lens
    effective_focal_length_m = -inv(c)
    back_focal_length_m = -a / c

    from_lyot = _multiply_abcd(
        lens,
        _translation_matrix(config.lyot_to_lens_m),
    )
    _, b_object, _, d_object = from_lyot
    conjugate_distance_m = -b_object / d_object
    conjugate = _multiply_abcd(
        _translation_matrix(conjugate_distance_m),
        from_lyot,
    )
    lateral_magnification = conjugate[1]
    detector_distance_m =
        config.lens_to_filter_m +
        config.filter_thickness_m / indices.fused_silica +
        config.filter_to_detector_m
    return (;
        indices,
        effective_focal_length_m,
        back_focal_length_m,
        conjugate_distance_m,
        detector_distance_m,
        detector_defocus_m=detector_distance_m - conjugate_distance_m,
        lateral_magnification,
    )
end

"""
    PreparedLLOWFS

Own the static optical assets, PROPER workspace, reflected Lyot coefficient,
LLOWFS relay aperture, scratch storage, and output intensity for repeated
propagation to the GoldEye sensor plane.
"""
struct PreparedLLOWFS{C,LC,E,CTX,W,AP,A,F,L,S,O,I} <: AbstractPreparedSpiders
    config::C
    llowfs_config::LC
    entrance::E
    context::CTX
    wavefront::W
    apertures::AP
    apodizer_internal::A
    fpm_opd_internal::F
    lyot_reflection_internal::L
    resample_scratch::S
    output::O
    refractive_indices::I
end

function _prepare_llowfs_relay!(
    wf::Proper.WaveFront,
    config::LLOWFSConfig,
    indices,
    context::Proper.RunContext,
)
    T = real(eltype(wf.field))
    prop_propagate(wf, T(config.lyot_to_lens_m), context, "LLOWFS lens front")
    front_focal_length = inv(
        T(indices.s_lah71 - 1) * T(config.front_curvature_inv_m),
    )
    prop_lens(wf, front_focal_length, context, "LLOWFS lens surface 1")
    lens = _prepare_clear_aperture_internal(
        wf,
        T(2 * config.lens_clear_radius_m),
    )
    _apply_internal_multiplier!(wf, lens)
    prop_propagate(
        wf,
        T(config.s_lah71_thickness_m / indices.s_lah71),
        context,
        "LLOWFS S-LAH71",
    )
    cement_focal_length = inv(
        T(indices.s_nph5 - indices.s_lah71) *
        T(config.cement_curvature_inv_m),
    )
    prop_lens(wf, cement_focal_length, context, "LLOWFS lens surface 2")
    prop_propagate(
        wf,
        T(_cement_to_detector_reduced_m(config, indices)),
        context,
        "LLOWFS lens surface 2 to GoldEye sensor",
    )
    return lens
end

function _propagate_llowfs_relay!(
    wf::Proper.WaveFront,
    config::LLOWFSConfig,
    indices,
    lens,
    context::Proper.RunContext,
)
    T = real(eltype(wf.field))
    prop_propagate(wf, T(config.lyot_to_lens_m), context, "LLOWFS lens front")
    front_focal_length = inv(
        T(indices.s_lah71 - 1) * T(config.front_curvature_inv_m),
    )
    _prepared_lens_aperture!(
        wf,
        front_focal_length,
        lens,
        context,
        "LLOWFS lens surface 1",
    )
    prop_propagate(
        wf,
        T(config.s_lah71_thickness_m / indices.s_lah71),
        context,
        "LLOWFS S-LAH71",
    )
    cement_focal_length = inv(
        T(indices.s_nph5 - indices.s_lah71) *
        T(config.cement_curvature_inv_m),
    )
    prop_lens(wf, cement_focal_length, context, "LLOWFS lens surface 2")
    prop_propagate(
        wf,
        T(_cement_to_detector_reduced_m(config, indices)),
        context,
        "LLOWFS lens surface 2 to GoldEye sensor",
    )
    return wf
end

"""
    prepare_llowfs(wavelength_m, gridsize; kwargs...)

Prepare the reflected-Lyot LLOWFS optical path for repeated propagation. The
entrance-array arguments have the same ownership and units as
[`prepare_spiders`](@ref). Detector pixel integration, noise, and hardware ROI
readout remain the responsibility of AdaptiveOpticsSim and the HIL model.
"""
function prepare_llowfs(
    wavelength_m::Real,
    gridsize::Integer;
    config::SpidersConfig=SpidersConfig(),
    llowfs_config::LLOWFSConfig=LLOWFSConfig(),
    ao_residual_nm::Union{Nothing,AbstractMatrix}=nothing,
    ao_sampling_m::Union{Nothing,Real}=nothing,
    pupil_amplitude::Union{Nothing,AbstractMatrix}=nothing,
    pupil_sampling_m::Union{Nothing,Real}=nothing,
    opd_m::Union{Nothing,AbstractMatrix}=nothing,
    opd_sampling_m::Union{Nothing,Real}=nothing,
    T::Type{<:AbstractFloat}=Float64,
)
    validate_prepared_config(config)
    validate(llowfs_config)
    wavelength_m > 0 || throw(ArgumentError("wavelength_m must be positive"))
    gridsize > 0 || throw(ArgumentError("gridsize must be positive"))
    T in (Float32, Float64) || throw(ArgumentError(
        "prepared propagation supports Float32 or Float64"))

    n = Int(gridsize)
    field = Matrix{Complex{T}}(undef, n, n)
    context = Proper.RunContext(typeof(field))
    wf = prop_begin!(
        field,
        config.telescope_diameter_m,
        T(wavelength_m);
        beam_diam_fraction=config.beam_diameter_fraction,
        context,
    )
    entrance = _prepared_entrance(
        wf,
        config,
        ao_residual_nm,
        ao_sampling_m,
        pupil_amplitude,
        pupil_sampling_m,
        opd_m,
        opd_sampling_m,
    )
    copyto!(wf.field, entrance.field)

    entrance_apertures = _prepare_to_apodizer!(wf, config, context)
    apodizer_internal = _prepare_apodizer_internal(wf, config, context)
    _apply_internal_multiplier!(wf, apodizer_internal)
    _propagate_to_fpm!(wf, context)
    fpm_opd_internal = _prepare_fpm_opd_internal(wf, config)
    _apply_internal_phase!(wf, fpm_opd_internal)
    lens0 = _prepare_to_lyot!(wf, context)
    lyot_reflection_internal = _prepare_lyot_internal(
        wf,
        config,
        Val(:reflection),
    )
    _apply_internal_multiplier!(wf, lyot_reflection_internal)
    indices = llowfs_refractive_indices(
        wavelength_m;
        temperature_c=llowfs_config.temperature_c,
    )
    llowfs_lens = _prepare_llowfs_relay!(wf, llowfs_config, indices, context)
    apertures = merge(entrance_apertures, (; lens0, llowfs_lens))

    resample_scratch = Matrix{T}(undef, n, n)
    output = Matrix{T}(undef, n, n)
    _prepared_intensity!(output, wf)
    prepared = PreparedLLOWFS(
        config,
        llowfs_config,
        entrance,
        context,
        wf,
        apertures,
        apodizer_internal,
        fpm_opd_internal,
        lyot_reflection_internal,
        resample_scratch,
        output,
        indices,
    )
    llowfs_propagate!(prepared)
    return prepared
end

"""
    llowfs_propagate!(prepared)

Propagate one frame through the reflected Lyot mask and Zemax-derived LLOWFS
relay. The centered sensor-plane intensity is written into and returned from
`prepared.output`. The warmed CPU call has a zero-byte Julia heap-allocation
contract; detector sampling and acquisition are outside this boundary.
"""
function llowfs_propagate!(prepared::PreparedLLOWFS)
    wf = _reset_prepared_entrance!(prepared, prepared.entrance)
    _propagate_to_apodizer!(
        wf,
        prepared.config,
        prepared.context,
        prepared.apertures,
    )
    _apply_internal_multiplier!(wf, prepared.apodizer_internal)
    _propagate_to_fpm!(wf, prepared.context)
    _apply_internal_phase!(wf, prepared.fpm_opd_internal)
    _propagate_to_lyot!(wf, prepared.context, prepared.apertures.lens0)
    _apply_internal_multiplier!(wf, prepared.lyot_reflection_internal)
    _propagate_llowfs_relay!(
        wf,
        prepared.llowfs_config,
        prepared.refractive_indices,
        prepared.apertures.llowfs_lens,
        prepared.context,
    )
    _prepared_intensity!(prepared.output, wf)
    return prepared.output
end

llowfs_intensity(prepared::PreparedLLOWFS) = prepared.output
