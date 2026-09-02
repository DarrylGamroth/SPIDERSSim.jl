using AdaptiveOpticsProperHIL
using AdaptiveOpticsSim
using FITSIO
using LinearAlgebra
using Proper
using SPIDERSSim
using SPIDERSSim.ProperOptics
using Test

using AdaptiveOpticsSim.AlgorithmGraphs

include("spiders_provisional.jl")
include("proper_optics.jl")
include("bax307.jl")
include("spiders_prescriptions.jl")
include("spiders_chopper.jl")
include("scc_pairing.jl")
include("reconstructors.jl")
