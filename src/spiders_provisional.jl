@enum SpidersEvidenceLevel::UInt8 begin
    SpidersDocumented = 1
    SpidersDeploymentConfigured = 2
    SpidersInferred = 3
    SpidersPlaceholder = 4
end

"""One reviewable claim behind a provisional SPIDERS instrument profile."""
struct SpidersProfileClaim
    name::Symbol
    evidence::SpidersEvidenceLevel
    current_choice::String
    source::String
    qualification_question::String
end

struct SpidersFocalPlaneMask{T<:Union{Float32,Float64}}
    label::Symbol
    stage_serial::String
    stage_counts::Int
    center_wavelength_um::T
    diameter_lambda_over_d::T
    vortex_charge::Int
end

struct SpidersGalilAlignment{T<:Union{Float32,Float64}}
    fpm_x_mm::T
    fpm_y_mm::T
    oae1_x_mm::T
    oae1_y_mm::T
    dm_x_mm::T
    dm_y_mm::T
end

struct SpidersFilterSelection{T<:Union{Float32,Float64}}
    wheel::Int
    serial::String
    slot::Int
    label::String
    center_wavelength_um::T
    full_width_um::T
end

struct SpidersOpticalGeometry{T<:Union{Float32,Float64}}
    subaru_pupil_diameter_m::T
    dm_pupil_diameter_m::T
    apodizer_pupil_diameter_m::T
    apodizer_outer_diameter_ratio::T
    apodizer_inner_diameter_ratio::T
    lyot_pupil_diameter_m::T
    lyot_outer_diameter_ratio::T
    lyot_inner_diameter_ratio::T
    lyot_spider_width_ratio::T
    focal_plane_f_number::T
    l2_focal_length_m::T
    fpm_to_l2_m::T
    reference_pinhole_diameter_m::T
    reference_pinhole_separation_pupils::T
    reference_pinhole_position_angle_deg::T
    camera_pixel_pitch_m::T
    camera_rotation_deg::T
end

"""
Cold, explicitly provisional SPIDERS H-band instrument profile.

Fields backed by the optical-design paper or SpiderMan deployment presets are
kept separate from the unresolved choices returned by
[`spiders_profile_claims`](@ref). This profile must not be treated as an
as-built or optically qualified prescription.
"""
struct ProvisionalSpidersProfile{T<:Union{Float32,Float64}}
    mode::Symbol
    wavelength_um::T
    focal_plane_mask::SpidersFocalPlaneMask{T}
    galil_alignment::SpidersGalilAlignment{T}
    llowfs_filter::SpidersFilterSelection{T}
    scc_filter::SpidersFilterSelection{T}
    llowfs_plane::Symbol
    llowfs_magnification::Symbol
    scc_plane::Symbol
    scc_main_beam::Symbol
    scc_reference_beam::Symbol
    llowfs_output_shape::Tuple{Int,Int}
    scc_output_shape::Tuple{Int,Int}
    optics::SpidersOpticalGeometry{T}
end

const _SPIDERS_WHEEL_SERIALS = (
    "TP02692146-23357",
    "TP02659529-22864",
)

function _spiders_h_filter(::Type{T}, wheel::Integer) where {
    T<:Union{Float32,Float64}
}
    wheel in (2, 3) || throw(ArgumentError(
        "a provisional SPIDERS science filter must use wheel 2 or 3",
    ))
    return SpidersFilterSelection{T}(
        Int(wheel),
        _SPIDERS_WHEEL_SERIALS[Int(wheel) - 1],
        5,
        "1550-25nm",
        T(1.550),
        T(0.025),
    )
end

function _spiders_provisional_optics(
    ::Type{T};
    reference_pinhole_diameter_m::Real,
    reference_pinhole_position_angle_deg::Real,
) where {T<:Union{Float32,Float64}}
    pinhole = T(reference_pinhole_diameter_m)
    isfinite(pinhole) && pinhole > zero(T) || throw(ArgumentError(
        "reference_pinhole_diameter_m must be finite and positive",
    ))
    angle = T(reference_pinhole_position_angle_deg)
    isfinite(angle) || throw(ArgumentError(
        "reference_pinhole_position_angle_deg must be finite",
    ))
    return SpidersOpticalGeometry{T}(
        T(7.920),
        T(0.033),
        T(0.012),
        T(7840.8 / 7920),
        T(2379.2 / 7920),
        T(0.00407),
        T(7527.2 / 7920),
        T(2645.7 / 7920),
        T(500 / 7920),
        T(64),
        T(0.286),
        T(0.354),
        pinhole,
        T(1.6),
        angle,
        T(15e-6),
        T(31),
    )
end

"""
    provisional_spiders_h_regular_profile(T=Float32; kwargs...)

Build the current review profile for the regular 1550 nm SPIDERS mask. The
default wheel ownership, selector states, SCC reference-hole choice, and
reference-hole position angle are placeholders. Every such choice is listed
by [`spiders_profile_claims`](@ref).
"""
function provisional_spiders_h_regular_profile(
    ::Type{T}=Float32;
    llowfs_wheel::Integer=2,
    scc_wheel::Integer=3,
    llowfs_plane::Symbol=:focal_plane,
    llowfs_magnification::Symbol=:zoomed_in,
    scc_plane::Symbol=:focal_plane,
    scc_main_beam::Symbol=:unblocked,
    scc_reference_beam::Symbol=:open,
    reference_pinhole_diameter_m::Real=203e-6,
    reference_pinhole_position_angle_deg::Real=0,
) where {T<:Union{Float32,Float64}}
    llowfs_wheel != scc_wheel || throw(ArgumentError(
        "the provisional LLOWFS and SCC profiles must use distinct wheels",
    ))
    llowfs_plane in (:pupil_plane, :focal_plane) || throw(ArgumentError(
        "llowfs_plane must be :pupil_plane or :focal_plane",
    ))
    llowfs_magnification in (:zoomed_out, :zoomed_in) || throw(ArgumentError(
        "llowfs_magnification must be :zoomed_out or :zoomed_in",
    ))
    scc_plane in (:pupil_plane, :focal_plane) || throw(ArgumentError(
        "scc_plane must be :pupil_plane or :focal_plane",
    ))
    scc_main_beam in (:blocked, :unblocked) || throw(ArgumentError(
        "scc_main_beam must be :blocked or :unblocked",
    ))
    scc_reference_beam in (:closed, :open) || throw(ArgumentError(
        "scc_reference_beam must be :closed or :open",
    ))

    return ProvisionalSpidersProfile{T}(
        :h_regular,
        T(1.550),
        SpidersFocalPlaneMask{T}(
            :h_regular,
            "27260961",
            424_960,
            T(1.550),
            T(5.7),
            4,
        ),
        SpidersGalilAlignment{T}(
            T(6.91),
            T(7.64),
            T(7.007),
            T(6.8435),
            T(6.054),
            T(6.827),
        ),
        _spiders_h_filter(T, llowfs_wheel),
        _spiders_h_filter(T, scc_wheel),
        llowfs_plane,
        llowfs_magnification,
        scc_plane,
        scc_main_beam,
        scc_reference_beam,
        (34, 47),
        (416, 380),
        _spiders_provisional_optics(
            T;
            reference_pinhole_diameter_m,
            reference_pinhole_position_angle_deg,
        ),
    )
end

"""
Return the unresolved or inferred choices that prevent optical qualification.

The tuple is deliberately cold metadata. Numerical propagation consumes the
typed profile fields directly and does not inspect these strings per frame.
"""
function spiders_profile_claims(profile::ProvisionalSpidersProfile)
    wheel_choice = "wheel $(profile.llowfs_filter.wheel) for LLOWFS; " *
        "wheel $(profile.scc_filter.wheel) for SCC"
    filter_choice = "slot 5 ($(profile.llowfs_filter.label)) on both arms"
    return (
        SpidersProfileClaim(
            :science_wheel_ownership,
            SpidersPlaceholder,
            wheel_choice,
            "SpiderMan deployment config lists two identical science-wheel " *
            "slot tables but does not identify their optical arms",
            "Which wheel serial is physically in the LLOWFS arm and which is " *
            "in the SCC arm?",
        ),
        SpidersProfileClaim(
            :science_filter_slot,
            SpidersInferred,
            filter_choice,
            "SpiderMan slot tables plus /mnt/datadrive/DATA/LOWFS/filt5",
            "Confirm the normal H-band slot and whether both arms use it.",
        ),
        SpidersProfileClaim(
            :llowfs_selector_state,
            SpidersPlaceholder,
            "$(profile.llowfs_plane), $(profile.llowfs_magnification)",
            "SpiderMan exposes both selector mechanisms but no qualified " *
            "observing-state preset",
            "Which LLOWFS plane and magnification are used for regular H?",
        ),
        SpidersProfileClaim(
            :scc_selector_state,
            SpidersPlaceholder,
            "$(profile.scc_plane), main $(profile.scc_main_beam), reference " *
            "$(profile.scc_reference_beam)",
            "SpiderMan exposes SCC pupil-lens and main-blocker mechanisms; " *
            "the acquisition service does not record one static state",
            "Confirm the SCC lens, main-beam, and chopper sequence for HIL.",
        ),
        SpidersProfileClaim(
            :reference_pinhole,
            SpidersPlaceholder,
            "$(profile.optics.reference_pinhole_diameter_m * 1e6) um at " *
            "$(profile.optics.reference_pinhole_separation_pupils) pupil " *
            "diameters and $(profile.optics.reference_pinhole_position_angle_deg) deg",
            "The optical design reports 203 and 291 um holes at 1.6 pupil " *
            "diameters, but the installed choice and position angle are absent",
            "Which reference hole is installed, and what is its pupil-plane " *
            "position angle?",
        ),
        SpidersProfileClaim(
            :focal_plane_mask_map,
            SpidersPlaceholder,
            "charge-4, 5.7 lambda/D regular-H geometry surrogate",
            "The optical design specifies mask diameter and family; no " *
            "as-built complex transmission map was found",
            "Provide the as-built tilt-Gaussian-vortex amplitude/phase map, " *
            "orientation, and registration.",
        ),
        SpidersProfileClaim(
            :apodizer_map,
            SpidersPlaceholder,
            "binary annular geometry surrogate",
            "The optical design gives normalized stops and a 10 um-dot " *
            "halftone description, but not the as-built transmission map",
            "Provide the as-built apodizer transmission and registration.",
        ),
        SpidersProfileClaim(
            :detector_mapping,
            SpidersDeploymentConfigured,
            "GoldEye LLOWFS $(profile.llowfs_output_shape) and C-RED 2 SCC " *
            "$(profile.scc_output_shape) products",
            "deployed image-adaptor streams, calibration scripts, and " *
            "/mnt/datadrive/DATA/Goldeye and CRED2-SCC products establish " *
            "camera ownership, array order, element type, and product shape",
            "Provide the qualified camera crops, optical sampling, " *
            "orientation, gains, distortion, and LLOWFS defocus/magnification.",
        ),
    )
end

spiders_profile_is_provisional(::ProvisionalSpidersProfile) = true
