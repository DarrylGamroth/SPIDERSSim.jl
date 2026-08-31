const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))

asset_path(name::AbstractString) = joinpath(PACKAGE_ROOT, name)

const DEFAULT_APODIZER_POLYNOMIAL = Float64[
    -4.6167316353135584e31,
    6.145914297871727e32,
    4.5600030735614572e32,
    -6.8907507359001675e32,
    -3.1347208721635807e32,
    -5.6520198759155703e32,
    -5.0805214834037565e32,
    6.0607863664877544e32,
    -2.0897545805393842e32,
    5.1067481423987977e32,
    -7.8025351459304506e31,
    6.9880517657408778e32,
    -1.9725110582963142e32,
    5.7591852922436367e32,
    -3.9085421357617009e32,
    -4.0351795617054639e32,
]

"""
    SpidersConfig(; kwargs...)

Configure the SPIDERS SCC optical prescription. Distances are in meters and
angles are in radians unless a field name states otherwise.

The defaults form a self-contained run using the bundled Subaru pupil and
apodizer FITS maps. Optional measured OAE maps and the deformable-mirror
correction are disabled because those calibration files are not bundled.
"""
Base.@kwdef struct SpidersConfig
    telescope_diameter_m::Float64 = 7.92
    input_fratio::Float64 = 13.901
    beam_diameter_fraction::Float64 = 0.1

    pupil_mode::Symbol = :fits
    pupil_path::String = asset_path("pupilsbr_nPup1200_kpdiam100_kodiam100_kthick100.fits")
    pupil_obscuration_ratio::Float64 = 0.29
    include_spiders::Bool = true

    turbulence::Bool = false
    ao_residual_path::String = asset_path("gpi2AoRes1arcsecSeeingMag8Snapshot.fits")
    scintillation::Bool = false

    oae1_error_path::Union{Nothing,String} = nothing
    oae2_error_path::Union{Nothing,String} = nothing
    dm_correction::Bool = false
    dm_actuators::Int = 24
    dm_actuators_across_pupil::Int = 23

    coronagraph::Bool = true
    apodizer_mode::Symbol = :fits
    apodizer_path::String = asset_path("pupilsbr_nPup1200_pdiam792_odiam230_thick000_Apod_rMask265.fits")
    apodizer_polynomial::Vector{Float64} = copy(DEFAULT_APODIZER_POLYNOMIAL)

    fpm_mode::Symbol = :analytic
    fpm_type::Symbol = :TG
    fpm_band::Symbol = :J
    fpm_map_path::Union{Nothing,String} = nothing
    fpm_map_sampling_m::Float64 = 250e-9
    fpm_map_type::Symbol = :mirror_surface
    center_pupils::Bool = true

    reference_pinhole::Bool = false
    lyot_pupil_diameter_m::Float64 = 2 * 2.034e-3
    lyot_pinhole_ratio::Float64 = 1 / 20
    lyot_reference_separation::Float64 = 1.6
    lyot_reference_angle_rad::Float64 = -pi / 2
    lyot_outer_margin::Float64 = 0.95
    lyot_obscuration_margin::Float64 = 1.15
    lyot_spider_margin::Float64 = 2.0
end

function validate(config::SpidersConfig)
    config.telescope_diameter_m > 0 || throw(ArgumentError("telescope_diameter_m must be positive"))
    config.input_fratio > 0 || throw(ArgumentError("input_fratio must be positive"))
    0 < config.beam_diameter_fraction <= 1 || throw(ArgumentError("beam_diameter_fraction must lie in (0, 1]"))
    0 <= config.pupil_obscuration_ratio < 1 || throw(ArgumentError("pupil_obscuration_ratio must lie in [0, 1)"))
    config.pupil_mode in (:fits, :analytic) || throw(ArgumentError("pupil_mode must be :fits or :analytic"))
    config.apodizer_mode in (:fits, :polynomial, :none) || throw(ArgumentError("apodizer_mode must be :fits, :polynomial, or :none"))
    config.fpm_mode in (:analytic, :map) || throw(ArgumentError("fpm_mode must be :analytic or :map"))
    config.fpm_type in (:TG, :TGV) || throw(ArgumentError("fpm_type must be :TG or :TGV"))
    config.fpm_band in (:J, :H) || throw(ArgumentError("the analytic FPM supports :J or :H"))
    config.fpm_map_type in (:phase, :wavefront, :mirror_surface) || throw(ArgumentError("invalid fpm_map_type"))
    config.dm_actuators > 0 || throw(ArgumentError("dm_actuators must be positive"))
    config.dm_actuators_across_pupil > 0 || throw(ArgumentError("dm_actuators_across_pupil must be positive"))
    config.lyot_pupil_diameter_m > 0 || throw(ArgumentError("lyot_pupil_diameter_m must be positive"))
    config.lyot_pinhole_ratio >= 0 || throw(ArgumentError("lyot_pinhole_ratio must be nonnegative"))

    config.pupil_mode === :fits && require_file(config.pupil_path, "pupil FITS map")
    config.apodizer_mode === :fits && require_file(config.apodizer_path, "apodizer FITS map")
    config.turbulence && require_file(config.ao_residual_path, "AO residual FITS map")
    config.oae1_error_path === nothing || require_file(config.oae1_error_path, "OAE1 error map")
    config.oae2_error_path === nothing || require_file(config.oae2_error_path, "OAE2 error map")
    if config.fpm_mode === :map
        config.fpm_map_path === nothing && throw(ArgumentError("fpm_map_path is required when fpm_mode=:map"))
        require_file(config.fpm_map_path, "FPM map")
    end
    return config
end

function require_file(path::AbstractString, description::AbstractString)
    isfile(path) || throw(ArgumentError("$description not found: $path"))
    return path
end
