abstract type AbstractSpidersEntrance end

struct StaticSpidersEntrance{F<:AbstractMatrix{<:Complex}} <: AbstractSpidersEntrance
    field::F
end

struct ResidualSpidersEntrance{
    F<:AbstractMatrix{<:Complex},
    R<:AbstractMatrix,
    T<:AbstractFloat,
} <: AbstractSpidersEntrance
    field::F
    residual_nm::R
    sampling_m::T
end

struct ExternalSpidersEntrance{
    F<:AbstractMatrix{<:Complex},
    P<:AbstractMatrix,
    O<:AbstractMatrix,
    T<:AbstractFloat,
} <: AbstractSpidersEntrance
    field::F
    pupil_amplitude::P
    opd_m::O
    pupil_sampling_m::T
    opd_sampling_m::T
end

"""
    PreparedSpiders

Own the static optical assets, PROPER workspace, wavefront, scratch storage,
and output intensity for repeated SPIDERS propagation. Construct it with
[`prepare_spiders`](@ref) and execute it with [`spiders_propagate!`](@ref).
"""
struct PreparedSpiders{C,E,CTX,W,AP,A,F,L,S,O}
    config::C
    entrance::E
    context::CTX
    wavefront::W
    apertures::AP
    apodizer_internal::A
    fpm_opd_internal::F
    lyot_internal::L
    resample_scratch::S
    output::O
end

function validate_prepared_config(config::SpidersConfig)
    validate(config)
    config.scintillation && throw(ArgumentError(
        "prepared propagation does not support scintillation yet"))
    config.turbulence && throw(ArgumentError(
        "bind an in-memory residual or OPD map instead of config.turbulence"))
    config.oae1_error_path === nothing || throw(ArgumentError(
        "prepared propagation does not support an OAE1 error map yet"))
    config.oae2_error_path === nothing || throw(ArgumentError(
        "prepared propagation does not support an OAE2 error map yet"))
    config.dm_correction && throw(ArgumentError(
        "prepared propagation expects the AO model to supply residual OPD"))
    config.fpm_mode === :analytic || throw(ArgumentError(
        "prepared propagation currently supports only fpm_mode=:analytic"))
    return config
end

function _base_entrance_field!(
    wf::Proper.WaveFront,
    config::SpidersConfig;
    external_pupil::Bool,
)
    prop_begin!(
        wf,
        config.telescope_diameter_m,
        prop_get_wavelength(wf);
        beam_diam_fraction=config.beam_diameter_fraction,
    )
    if external_pupil
        prop_circular_aperture(wf, config.telescope_diameter_m / 2)
        prop_define_entrance(wf)
    else
        apply_entrance_pupil!(wf, config)
    end
    return copy(wf.field)
end

function _prepared_entrance(
    wf::Proper.WaveFront,
    config::SpidersConfig,
    ao_residual_nm,
    ao_sampling_m,
    pupil_amplitude,
    pupil_sampling_m,
    opd_m,
    opd_sampling_m,
)
    has_external_pupil = pupil_amplitude !== nothing
    has_external_opd = opd_m !== nothing
    has_external_pupil == has_external_opd || throw(ArgumentError(
        "pupil_amplitude and opd_m must be supplied together"))
    ao_residual_nm === nothing || !has_external_opd || throw(ArgumentError(
        "ao_residual_nm cannot be combined with pupil_amplitude/opd_m"))

    field = _base_entrance_field!(wf, config;
        external_pupil=has_external_pupil)
    T = real(eltype(wf.field))
    if has_external_pupil
        pupil_sampling_m === nothing && throw(ArgumentError(
            "pupil_sampling_m is required with pupil_amplitude"))
        resolved_opd_sampling = something(opd_sampling_m, pupil_sampling_m)
        return ExternalSpidersEntrance(
            field,
            pupil_amplitude,
            opd_m,
            T(pupil_sampling_m),
            T(resolved_opd_sampling),
        )
    elseif ao_residual_nm !== nothing
        sampling = something(
            ao_sampling_m,
            0.045 / 7.73 * config.telescope_diameter_m,
        )
        return ResidualSpidersEntrance(field, ao_residual_nm, T(sampling))
    end
    return StaticSpidersEntrance(field)
end

@inline function _apply_internal_multiplier!(
    wf::Proper.WaveFront,
    multiplier::AbstractMatrix,
)
    @. wf.field *= multiplier
    return wf
end

@inline _apply_internal_multiplier!(wf::Proper.WaveFront, ::Nothing) = wf

@inline function _apply_internal_phase!(
    wf::Proper.WaveFront,
    opd_m::AbstractMatrix,
)
    scale = 2pi / wf.wavelength_m
    @. wf.field *= cis(scale * opd_m)
    return wf
end

@inline _apply_internal_phase!(wf::Proper.WaveFront, ::Nothing) = wf

function _centered_to_internal(wf::Proper.WaveFront, values::AbstractMatrix)
    T = real(eltype(wf.field))
    internal = similar(wf.field, T, size(wf.field))
    prop_shift_center!(internal, values; inverse=true)
    return internal
end

@inline function _prepared_intensity!(
    output::Matrix,
    wf::Proper.WaveFront{T,F},
) where {T,F<:Matrix}
    field = wf.field
    ny, nx = size(field)
    size(output) == (ny, nx) || throw(DimensionMismatch(
        "prepared intensity output must match the wavefront grid"))
    sy = ny ÷ 2
    sx = nx ÷ 2
    @inbounds for column in 1:nx
        source_column = mod1(column + sx, nx)
        for row in 1:ny
            source_row = mod1(row + sy, ny)
            output[row, column] = abs2(field[source_row, source_column])
        end
    end
    return output
end

@inline _prepared_intensity!(output::AbstractMatrix, wf::Proper.WaveFront) =
    prop_end!(output, wf)

function _prepare_apodizer_internal(
    wf::Proper.WaveFront,
    config::SpidersConfig,
    context::Proper.RunContext,
)
    (!config.coronagraph || config.apodizer_mode === :none) && return nothing
    if config.apodizer_mode === :fits
        source = prop_readmap(
            wf,
            config.apodizer_path,
            context;
            SAMPLING=2 * prop_get_beamradius(wf) / 1200,
        )
        internal = similar(wf.field, real(eltype(wf.field)), size(wf.field))
        copyto!(internal, source)
        return internal
    end
    centered = radial_apodizer(
        config.apodizer_polynomial,
        prop_get_gridsize(wf),
        prop_get_beamradius(wf),
        config.pupil_obscuration_ratio,
        prop_get_sampling(wf),
    )
    return _centered_to_internal(wf, centered)
end

function _prepare_fpm_opd_internal(
    wf::Proper.WaveFront,
    config::SpidersConfig,
)
    config.coronagraph || return nothing
    wavelength = prop_get_wavelength(wf)
    radius = config.fpm_band === :J ? 459.62e-6 / 2 : 584.63e-6 / 2
    phase = if config.fpm_type === :TG
        tilt_gaussian_fpm(
            wf;
            radius_m=radius,
            tilt_rad=1.6 * atan(1 / 64),
            gaussian_phase=7.8,
            charge=0,
            piston_rad=1.9999pi,
            levels=16,
            center_pupils=config.center_pupils,
        )
    else
        sampling = prop_get_sampling(wf)
        angular_sampling = prop_get_sampling_radians(wf)
        resolution_element_m = wavelength / wf.beam_diameter_m /
                               angular_sampling * sampling
        tilt_gaussian_fpm(
            wf;
            radius_m=resolution_element_m,
            tilt_rad=1.6 * atan(1 / 64),
            gaussian_phase=10,
            charge=4,
            piston_rad=0,
            levels=16,
            center_pupils=config.center_pupils,
        )
    end
    phase .*= wavelength / 2pi
    return _centered_to_internal(wf, phase)
end

function _prepare_lyot_internal(
    wf::Proper.WaveFront,
    config::SpidersConfig,
)
    config.coronagraph || return nothing
    pupil_diameter = config.lyot_pupil_diameter_m
    center_y = config.center_pupils ?
               config.lyot_reference_separation * pupil_diameter / 2 : 0.0
    pinhole_diameter = config.reference_pinhole ?
                       config.lyot_pinhole_ratio * pupil_diameter : 0.0
    lyot = scc_lyot_stop(
        wf,
        pupil_diameter;
        pinhole_diameter_m=pinhole_diameter,
        separation_m=config.lyot_reference_separation * pupil_diameter,
        position_angle_rad=config.lyot_reference_angle_rad,
        obscuration_diameter_m=config.pupil_obscuration_ratio * pupil_diameter,
        center=(0.0, center_y),
        outer_margin=config.lyot_outer_margin,
        obscuration_margin=config.lyot_obscuration_margin,
        spider_margin=config.lyot_spider_margin,
    )
    return _centered_to_internal(wf, lyot.mask)
end

@inline function _prepared_lens_aperture!(
    wf::Proper.WaveFront,
    focal_length_m::Real,
    aperture_internal::AbstractMatrix,
    context::Proper.RunContext,
    name::AbstractString,
)
    T = real(eltype(wf.field))
    prop_lens(wf, T(focal_length_m), context, name)
    _apply_internal_multiplier!(wf, aperture_internal)
    return wf
end

function _prepare_clear_aperture_internal(
    wf::Proper.WaveFront,
    diameter_m::Real,
)
    T = real(eltype(wf.field))
    radius = T(diameter_m) / T(2)
    centered = prop_ellipse(wf, radius, radius)
    return _centered_to_internal(wf, centered)
end

function _prepare_to_apodizer!(
    wf::Proper.WaveFront,
    config::SpidersConfig,
    context::Proper.RunContext,
)
    T = real(eltype(wf.field))
    diameter = T(config.telescope_diameter_m)
    input_focal_length = diameter * T(config.input_fratio)
    prop_propagate(wf, input_focal_length, context, "dummy lens")
    prop_lens(wf, input_focal_length, context, "F/13.901 convergent beam")
    prop_propagate(wf, input_focal_length, context, "Bay #4 focus")
    prop_propagate(wf, T(772.35e-3), context, "OAE1")
    prop_lens(wf, T(459.7e-3), context, "OAE1")
    oae1 = _prepare_clear_aperture_internal(wf, T(3 * 25.4e-3))
    _apply_internal_multiplier!(wf, oae1)
    prop_propagate(wf, T(471.4e-3), context, "DM")
    prop_propagate(
        wf,
        prop_get_distancetofocus(wf) + T(367.2e-3),
        context,
        "OAE2",
    )
    prop_lens(wf, T(277.8e-3), context, "OAE2")
    oae2 = _prepare_clear_aperture_internal(wf, T(38.1e-3))
    _apply_internal_multiplier!(wf, oae2)
    prop_propagate(wf, T(377.2e-3), context, "Apodizer")
    return (; oae1, oae2)
end

function _propagate_to_apodizer!(
    wf::Proper.WaveFront,
    config::SpidersConfig,
    context::Proper.RunContext,
    apertures,
)
    T = real(eltype(wf.field))
    diameter = T(config.telescope_diameter_m)
    input_focal_length = diameter * T(config.input_fratio)
    prop_propagate(wf, input_focal_length, context, "dummy lens")
    prop_lens(wf, input_focal_length, context, "F/13.901 convergent beam")
    prop_propagate(wf, input_focal_length, context, "Bay #4 focus")
    prop_propagate(wf, T(772.35e-3), context, "OAE1")
    _prepared_lens_aperture!(wf, T(459.7e-3), apertures.oae1, context, "OAE1")
    prop_propagate(wf, T(471.4e-3), context, "DM")
    prop_propagate(
        wf,
        prop_get_distancetofocus(wf) + T(367.2e-3),
        context,
        "OAE2",
    )
    _prepared_lens_aperture!(wf, T(277.8e-3), apertures.oae2, context, "OAE2")
    prop_propagate(wf, T(377.2e-3), context, "Apodizer")
    return wf
end

@inline function _propagate_to_fpm!(
    wf::Proper.WaveFront,
    context::Proper.RunContext,
)
    prop_propagate(wf, prop_get_distancetofocus(wf), context, "FPM")
    return wf
end

function _prepare_to_lyot!(
    wf::Proper.WaveFront,
    context::Proper.RunContext,
)
    T = real(eltype(wf.field))
    prop_propagate(wf, T(358.69e-3 + 5.166e-3), context, "Lens 0")
    prop_lens(wf, T(285.99e-3), context, "Lens 0")
    lens0 = _prepare_clear_aperture_internal(wf, T(40e-3))
    _apply_internal_multiplier!(wf, lens0)
    prop_propagate(wf, T(382.67e-3 - 0.022e-3), context, "Lyot")
    return lens0
end

function _propagate_to_lyot!(
    wf::Proper.WaveFront,
    context::Proper.RunContext,
    lens0,
)
    T = real(eltype(wf.field))
    prop_propagate(wf, T(358.69e-3 + 5.166e-3), context, "Lens 0")
    _prepared_lens_aperture!(wf, T(285.99e-3), lens0, context, "Lens 0")
    prop_propagate(wf, T(382.67e-3 - 0.022e-3), context, "Lyot")
    return wf
end

function _prepare_to_scc_focus!(
    wf::Proper.WaveFront,
    context::Proper.RunContext,
)
    T = real(eltype(wf.field))
    prop_propagate(wf, T(55e-3), context, "Chopper")
    prop_propagate(wf, T(133.16e-3 - 0.022e-3), context, "Lens 1")
    prop_lens(wf, T(285.99e-3), context, "Lens 1")
    lens1 = _prepare_clear_aperture_internal(wf, T(40e-3))
    _apply_internal_multiplier!(wf, lens1)
    prop_propagate(
        wf,
        prop_get_distancetofocus(wf),
        context,
        "SCC Focus",
    )
    return lens1
end


function _propagate_to_scc_focus!(
    wf::Proper.WaveFront,
    context::Proper.RunContext,
    lens1,
)
    T = real(eltype(wf.field))
    prop_propagate(wf, T(55e-3), context, "Chopper")
    prop_propagate(wf, T(133.16e-3 - 0.022e-3), context, "Lens 1")
    _prepared_lens_aperture!(wf, T(285.99e-3), lens1, context, "Lens 1")
    prop_propagate(
        wf,
        prop_get_distancetofocus(wf),
        context,
        "SCC Focus",
    )
    return wf
end

"""
    prepare_spiders(wavelength_m, gridsize; kwargs...)

Prepare the CPU SPIDERS optical path for repeated propagation. Static FITS
maps and focal/Lyot masks are loaded or constructed once. The returned owner
binds its optional input arrays by identity; update their contents between
calls rather than replacing them.

Use `ao_residual_nm` for MATLAB-compatible standalone runs. For integration
with AdaptiveOpticsSim, bind `pupil_amplitude`, `pupil_sampling_m`, `opd_m`,
and optionally `opd_sampling_m`; the OPD is in meters and replaces the bundled
entrance pupil and internal DM correction.
"""
function prepare_spiders(
    wavelength_m::Real,
    gridsize::Integer;
    config::SpidersConfig=SpidersConfig(),
    ao_residual_nm::Union{Nothing,AbstractMatrix}=nothing,
    ao_sampling_m::Union{Nothing,Real}=nothing,
    pupil_amplitude::Union{Nothing,AbstractMatrix}=nothing,
    pupil_sampling_m::Union{Nothing,Real}=nothing,
    opd_m::Union{Nothing,AbstractMatrix}=nothing,
    opd_sampling_m::Union{Nothing,Real}=nothing,
    T::Type{<:AbstractFloat}=Float64,
)
    validate_prepared_config(config)
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
    lyot_internal = _prepare_lyot_internal(wf, config)
    _apply_internal_multiplier!(wf, lyot_internal)
    lens1 = _prepare_to_scc_focus!(wf, context)
    apertures = merge(entrance_apertures, (; lens0, lens1))

    resample_scratch = Matrix{T}(undef, n, n)
    output = Matrix{T}(undef, n, n)
    _prepared_intensity!(output, wf)
    prepared = PreparedSpiders(
        config,
        entrance,
        context,
        wf,
        apertures,
        apodizer_internal,
        fpm_opd_internal,
        lyot_internal,
        resample_scratch,
        output,
    )
    spiders_propagate!(prepared)
    return prepared
end

@inline function _reset_prepared_field!(
    prepared::PreparedSpiders,
    field::AbstractMatrix{<:Complex},
)
    prop_begin!(
        prepared.wavefront,
        prepared.config.telescope_diameter_m,
        prop_get_wavelength(prepared.wavefront);
        beam_diam_fraction=prepared.config.beam_diameter_fraction,
    )
    copyto!(prepared.wavefront.field, field)
    return prepared.wavefront
end


@inline _reset_prepared_entrance!(
    prepared::PreparedSpiders,
    entrance::StaticSpidersEntrance,
) = _reset_prepared_field!(prepared, entrance.field)

@inline function _reset_prepared_entrance!(
    prepared::PreparedSpiders,
    entrance::ResidualSpidersEntrance,
)
    wf = _reset_prepared_field!(prepared, entrance.field)
    residual = entrance.residual_nm
    prop_resamplemap!(
        prepared.resample_scratch,
        wf,
        residual,
        entrance.sampling_m,
        size(residual, 2) / 2,
        size(residual, 1) / 2,
        prepared.context,
    )
    prepared.resample_scratch .*= 1e-9
    prop_add_phase(wf, prepared.resample_scratch)
    return wf
end

@inline function _reset_prepared_entrance!(
    prepared::PreparedSpiders,
    entrance::ExternalSpidersEntrance,
)
    wf = _reset_prepared_field!(prepared, entrance.field)
    pupil = entrance.pupil_amplitude
    prop_resamplemap!(
        prepared.resample_scratch,
        wf,
        pupil,
        entrance.pupil_sampling_m,
        size(pupil, 2) / 2,
        size(pupil, 1) / 2,
        prepared.context,
    )
    prop_multiply(wf, prepared.resample_scratch)

    opd = entrance.opd_m
    prop_resamplemap!(
        prepared.resample_scratch,
        wf,
        opd,
        entrance.opd_sampling_m,
        size(opd, 2) / 2,
        size(opd, 1) / 2,
        prepared.context,
    )
    prop_add_phase(wf, prepared.resample_scratch)
    return wf
end

"""
    spiders_propagate!(prepared)

Propagate one frame using the arrays bound by [`prepare_spiders`](@ref), write
the centered final intensity into `prepared.output`, and return that same
buffer. Compilation, preparation, input-array replacement, diagnostics, and
detector effects are outside this steady-state operation boundary.
"""
function spiders_propagate!(prepared::PreparedSpiders)
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
    _apply_internal_multiplier!(wf, prepared.lyot_internal)
    _propagate_to_scc_focus!(
        wf,
        prepared.context,
        prepared.apertures.lens1,
    )
    _prepared_intensity!(prepared.output, wf)
    return prepared.output
end

spiders_intensity(prepared::PreparedSpiders) = prepared.output
