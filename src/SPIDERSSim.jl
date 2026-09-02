"""Instrument-level SPIDERS optical simulation and RTC boundary."""
module SPIDERSSim

using AdaptiveOpticsProperHIL
using AdaptiveOpticsSim
using LinearAlgebra
using Proper

import AdaptiveOpticsProperHIL:
    prepare_proper_assets,
    prepare_proper_state,
    reset_proper_state!

const _AOG = AdaptiveOpticsSim.AlgorithmGraphs
const _Backends = AdaptiveOpticsSim.Backends

@inline _prop_lens!(wavefront, focal_length_m, ::Nothing, name) =
    Proper.prop_lens(wavefront, focal_length_m, name)
@inline _prop_lens!(
    wavefront,
    focal_length_m,
    run_context::Proper.RunContext,
    name,
) = Proper.prop_lens(wavefront, focal_length_m, run_context, name)

@inline _prop_propagate!(wavefront, distance_m, ::Nothing, name) =
    Proper.prop_propagate(wavefront, distance_m, name)
@inline _prop_propagate!(
    wavefront,
    distance_m,
    run_context::Proper.RunContext,
    name,
) = Proper.prop_propagate(wavefront, distance_m, run_context, name)

export SpidersEvidenceLevel,
    SpidersDocumented,
    SpidersDeploymentConfigured,
    SpidersInferred,
    SpidersPlaceholder,
    SpidersProfileClaim,
    SpidersFocalPlaneMask,
    SpidersGalilAlignment,
    SpidersFilterSelection,
    SpidersOpticalGeometry,
    ProvisionalSpidersProfile,
    provisional_spiders_h_regular_profile,
    provisional_spiders_entrance_pupil_amplitude,
    spiders_profile_claims,
    spiders_profile_is_provisional,
    SpidersSurrogateParameters,
    SpidersChopperPhase,
    SpidersUnfringed,
    SpidersFringed,
    SpidersChopperPlan,
    SpidersChopperFrameStatus,
    SpidersSCCPairPlan,
    SpidersSCCPairState,
    SpidersSCCPairProducts,
    prepare_scc_pairing,
    accept_scc_frame!,
    reset_scc_pairing!,
    ProvisionalSpidersLLOWFSPrescription,
    ProvisionalSpidersSCCPrescription,
    ProvisionalSpidersChoppedSCCPrescription,
    PROVISIONAL_SPIDERS_LLOWFS_SCHEMA,
    PROVISIONAL_SPIDERS_SCC_SCHEMA,
    PROVISIONAL_SPIDERS_CHOPPED_SCC_SCHEMA,
    provisional_spiders_llowfs_configuration,
    provisional_spiders_scc_configuration,
    provisional_spiders_chopped_scc_configuration,
    spiders_chopper_frame_status,
    spiders_chopper_sequence,
    spiders_chopper_phase,
    spiders_prescription_claims,
    spiders_configuration_claims,
    SpidersLLOWFSPlan,
    SpidersLLOWFSWorkspace,
    SpidersLLOWFSProducts,
    prepare_llowfs,
    process_llowfs!,
    SpidersSCCPlan,
    SpidersSCCWorkspace,
    SpidersSCCProducts,
    prepare_scc,
    process_scc!

export BAX307CalibrationPlan,
    BAX307_ACTUATOR_COUNT,
    BAX307_GRID_SIZE,
    BAX307_PDM_COMMAND_SCHEMA,
    BAX307_SURFACE_OPD_SCHEMA,
    BAX307_ACTUATOR_COORDINATES_SCHEMA,
    BAX307_INFLUENCE_FUNCTIONS_SCHEMA,
    prepare_bax307_calibration,
    load_bax307_calibration,
    load_bax307_best_flat,
    bax307_actuator_coordinates,
    bax307_influence_functions_opd,
    bax307_surface_opd,
    bax307_valid_command_mask,
    bax307_actuator_model,
    bax307_deformable_mirror_node,
    bax307_graph_parameters

include("spiders_provisional.jl")
include("bax307.jl")
include("spiders_prescriptions.jl")
include("spiders_chopper.jl")
include("scc_pairing.jl")
include("reconstructors.jl")

end # module SPIDERSSim
