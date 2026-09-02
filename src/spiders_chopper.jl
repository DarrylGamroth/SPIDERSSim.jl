const PROVISIONAL_SPIDERS_CHOPPED_SCC_SCHEMA =
    "org.subaru.spiders.scc-chopped-relative-intensity/provisional-1"

@enum SpidersChopperPhase::UInt8 begin
    SpidersUnfringed = 0
    SpidersFringed = 1
end

"""
Run-immutable contract for ideal camera-synchronous SCC chopping.

Each successful graph step is one complete exposure. The reference path is
fully open for a fringed exposure and fully closed for an unfringed exposure;
blade transitions are assumed to occur outside the modeled exposure. Model
time and camera cadence remain owned by the graph driver.
"""
struct SpidersChopperPlan
    initial_phase::SpidersChopperPhase
end

SpidersChopperPlan(;
    initial_phase::SpidersChopperPhase=SpidersFringed,
) = SpidersChopperPlan(initial_phase)

"""Identity metadata for one successfully completed chopped SCC frame."""
struct SpidersChopperFrameStatus
    sequence::UInt64
    phase::SpidersChopperPhase
end

@inline spiders_chopper_sequence(status::SpidersChopperFrameStatus) =
    status.sequence
@inline spiders_chopper_phase(status::SpidersChopperFrameStatus) = status.phase

struct ProvisionalSpidersChoppedSCCPrescription{
    T<:Union{Float32,Float64},
} <: Function
    profile::ProvisionalSpidersProfile{T}
    parameters::SpidersSurrogateParameters{T}
    chopper::SpidersChopperPlan
end

mutable struct _SpidersChopperState
    completed_sequence::UInt64
end

"""
    provisional_spiders_chopped_scc_configuration(profile; kwargs...)

Create the ideal synchronized-chopper variant of the provisional SCC Proper
configuration. One graph step produces one complete fringed or unfringed
frame. The initial phase defaults to `SpidersFringed`; every subsequent
successful step alternates phase. Use [`spiders_chopper_frame_status`](@ref)
after each step to attach the exact sequence and phase at the HIL boundary.
"""
function provisional_spiders_chopped_scc_configuration(
    profile::ProvisionalSpidersProfile{T};
    resolution::Integer=512,
    pupil_resolution::Union{Nothing,Integer}=nothing,
    beam_diameter_fraction::Real=0.25,
    camera_magnification::Real=1,
    detector_defocus_m::Real=0,
    central_tilt_cycles::Real=1,
    central_gaussian_phase_rad::Real=pi,
    central_gaussian_sigma_ratio::Real=0.35,
    initial_phase::SpidersChopperPhase=SpidersFringed,
    rng_seed::Integer=0,
    output_schema::AbstractString=PROVISIONAL_SPIDERS_CHOPPED_SCC_SCHEMA,
) where {T}
    profile.scc_reference_beam === :open || throw(ArgumentError(
        "a chopped SCC configuration requires scc_reference_beam=:open; " *
        "the chopper owns per-exposure reference transmission",
    ))
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
    prescription = ProvisionalSpidersChoppedSCCPrescription(
        profile,
        parameters,
        SpidersChopperPlan(; initial_phase),
    )
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

function spiders_prescription_claims(
    configuration::ProperPropagationConfiguration{
        T,
        <:ProvisionalSpidersChoppedSCCPrescription,
    },
) where {T}
    prescription = configuration.prescription
    detector_choice = "chopped SCC $(configuration.output_rows)-by-" *
        "$(configuration.output_columns), $(prescription.profile.scc_plane), " *
        "magnification $(prescription.parameters.camera_magnification), " *
        "defocus $(prescription.parameters.detector_defocus_m) m, rotation " *
        "$(prescription.profile.optics.camera_rotation_deg) deg"
    return (
        _spiders_prescription_claims(
            configuration,
            prescription,
            detector_choice,
        )...,
        SpidersProfileClaim(
            :frame_synchronous_chopper,
            SpidersPlaceholder,
            "one complete fringed exposure followed by one complete " *
            "unfringed exposure; initial phase " *
            "$(prescription.chopper.initial_phase)",
            "The simplechopper firmware phase-locks a ten-slot wheel to " *
            "camera-trigger and blade signals, but the exact exposure, duty, " *
            "phase, and transition margin are not qualified",
            "Confirm the camera trigger rate, exposure duration, initial " *
            "parity, wheel duty cycle, phase offset, and blade transition " *
            "margin.",
        ),
    )
end

function spiders_configuration_claims(
    configuration::ProperPropagationConfiguration{
        T,
        <:ProvisionalSpidersChoppedSCCPrescription,
    },
) where {T}
    return (
        spiders_profile_claims(configuration.prescription.profile)...,
        spiders_prescription_claims(configuration)...,
    )
end

@inline _spiders_output_shape(
    prescription::ProvisionalSpidersChoppedSCCPrescription,
) = prescription.profile.scc_output_shape

function _spiders_chopped_lyot_masks(
    prescription::ProvisionalSpidersChoppedSCCPrescription{T},
    main_lyot::Matrix{T},
    reference_hole::Matrix{T},
) where {T}
    main_scale = prescription.profile.scc_main_beam === :unblocked ?
        one(T) : zero(T)
    unfringed = similar(main_lyot)
    fringed = similar(main_lyot)
    @. unfringed = main_scale * main_lyot
    @. fringed = unfringed + reference_hole
    return fringed, unfringed
end

function prepare_proper_assets(
    prescription::ProvisionalSpidersChoppedSCCPrescription{T},
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
    fringed, unfringed = _spiders_chopped_lyot_masks(
        prescription,
        masks.main_lyot,
        masks.reference_hole,
    )
    return merge(assets, (
        fringed_lyot_mask=_spiders_target_array(target, fringed),
        unfringed_lyot_mask=_spiders_target_array(target, unfringed),
    ))
end

prepare_proper_state(
    ::ProvisionalSpidersChoppedSCCPrescription,
    ::Type,
) = _SpidersChopperState(UInt64(0))

function reset_proper_state!(
    ::ProvisionalSpidersChoppedSCCPrescription,
    state::_SpidersChopperState,
)
    state.completed_sequence = UInt64(0)
    return state
end

@inline _spiders_rotate_camera(
    prescription::ProvisionalSpidersChoppedSCCPrescription,
    square_intensity,
    rotated_intensity,
    run_context,
) = _spiders_rotate_scc_camera(
    prescription.profile,
    square_intensity,
    rotated_intensity,
    run_context,
)

@inline _spiders_camera_plane!(
    prescription::ProvisionalSpidersChoppedSCCPrescription,
    wavefront,
    run_context,
) = _spiders_scc_camera_plane!(
    prescription.profile,
    prescription.parameters,
    wavefront,
    run_context,
)

@inline function _spiders_other_chopper_phase(phase::SpidersChopperPhase)
    return phase === SpidersFringed ? SpidersUnfringed : SpidersFringed
end

@inline function _spiders_chopper_phase_at(
    plan::SpidersChopperPlan,
    zero_based_sequence::UInt64,
)
    return iseven(zero_based_sequence) ? plan.initial_phase :
        _spiders_other_chopper_phase(plan.initial_phase)
end

@noinline function _throw_spiders_chopper_sequence_exhausted()
    throw(ArgumentError("the SPIDERS chopper frame sequence is exhausted"))
end

@inline function _spiders_next_chopper_phase(
    prescription::ProvisionalSpidersChoppedSCCPrescription,
    state::_SpidersChopperState,
)
    state.completed_sequence == typemax(UInt64) &&
        _throw_spiders_chopper_sequence_exhausted()
    return _spiders_chopper_phase_at(
        prescription.chopper,
        state.completed_sequence,
    )
end

@inline function _spiders_commit_chopper_frame!(state::_SpidersChopperState)
    state.completed_sequence += UInt64(1)
    return state
end

@inline function _spiders_chopper_lyot_mask(
    phase::SpidersChopperPhase,
    fringed_lyot_mask,
    unfringed_lyot_mask,
)
    return phase === SpidersFringed ?
        fringed_lyot_mask : unfringed_lyot_mask
end

function (prescription::ProvisionalSpidersChoppedSCCPrescription)(
    wavelength_m,
    resolution;
    prescription_state,
    fringed_lyot_mask,
    unfringed_lyot_mask,
    kwargs...,
)
    phase = _spiders_next_chopper_phase(prescription, prescription_state)
    lyot_mask = _spiders_chopper_lyot_mask(
        phase,
        fringed_lyot_mask,
        unfringed_lyot_mask,
    )
    result = _run_provisional_spiders_prescription(
        prescription,
        wavelength_m,
        resolution;
        lyot_mask,
        kwargs...,
    )
    _spiders_commit_chopper_frame!(prescription_state)
    return result
end

@noinline function _throw_spiders_chopper_frame_unavailable()
    throw(ArgumentError(
        "no chopped SCC frame has completed since preparation or reset",
    ))
end

@inline function _spiders_chopper_frame_status(
    prescription::ProvisionalSpidersChoppedSCCPrescription,
    state::_SpidersChopperState,
)
    sequence = state.completed_sequence
    iszero(sequence) && _throw_spiders_chopper_frame_unavailable()
    phase = _spiders_chopper_phase_at(
        prescription.chopper,
        sequence - UInt64(1),
    )
    return SpidersChopperFrameStatus(sequence, phase)
end

function _spiders_chopper_frame_status(prescription, state)
    throw(ArgumentError(
        "the selected graph node is not a chopped SPIDERS SCC prescription",
    ))
end

@inline function _validate_spiders_chopper_graph_status(graph, status)
    _AOG.graph_failed(graph) && throw(ArgumentError(
        "a failed graph has no publishable chopped SCC frame",
    ))
    status.sequence == _AOG.graph_step_sequence(graph) || throw(ArgumentError(
        "the chopped SCC sequence is not aligned with the graph sequence",
    ))
    return status
end

@inline function spiders_chopper_frame_status(
    graph::_AOG.PreparedAlgorithmGraph,
    ::Val{Name},
) where {Name}
    owner = _AOG.prepared_graph_node(graph, Val(Name))
    status = _spiders_chopper_frame_status(
        owner.plan.prescription,
        owner.state.prescription,
    )
    return _validate_spiders_chopper_graph_status(graph, status)
end

@inline spiders_chopper_frame_status(graph::_AOG.PreparedAlgorithmGraph) =
    spiders_chopper_frame_status(graph, Val(:proper))

function spiders_chopper_frame_status(
    graph::_AOG.PreparedAlgorithmGraph,
    name::Symbol,
)
    owner = _AOG.prepared_graph_node(graph, name)
    status = _spiders_chopper_frame_status(
        owner.plan.prescription,
        owner.state.prescription,
    )
    return _validate_spiders_chopper_graph_status(graph, status)
end
