"""Instrument-level SPIDERS optical simulation and RTC boundary."""
module SPIDERSSim

using AdaptiveOpticsProperHIL
using AdaptiveOpticsSim
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
    spiders_profile_claims,
    spiders_profile_is_provisional,
    SpidersSurrogateParameters,
    SpidersChopperPhase,
    SpidersUnfringed,
    SpidersFringed,
    SpidersChopperPlan,
    SpidersChopperFrameStatus,
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
    spiders_configuration_claims

include("spiders_provisional.jl")
include("spiders_prescriptions.jl")
include("spiders_chopper.jl")

end # module SPIDERSSim

