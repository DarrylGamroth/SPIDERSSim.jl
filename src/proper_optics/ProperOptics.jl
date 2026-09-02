"""
SPIDERS SCC and reflected-Lyot LLOWFS physical-optics propagation.

The prescription was incorporated from SpidersProper revision
`09d92f5fd74d4b77869c7f430a06b3e077795691`. Reference optical products remain
external inputs rather than package source assets.
"""
module ProperOptics

using FFTW: fft, fftshift, ifft, ifftshift
using Proper
using Random: AbstractRNG, default_rng, rand, randn
using Statistics: median, quantile

export SCCPropagationConfiguration,
    SCCPropagationResult,
    PreparedSCCPropagation,
    prepare_scc_propagation,
    propagate_scc!,
    propagate_scc,
    propagate_scc_movie,
    scc_intensity,
    scc_field

export LLOWFSRelayConfiguration,
    PreparedLLOWFSPropagation,
    prepare_llowfs_propagation,
    propagate_llowfs!,
    llowfs_intensity,
    llowfs_refractive_indices,
    llowfs_relay_summary

export LyotStop,
    subaru_pupil!,
    radial_apodizer,
    tilt_gaussian_fpm,
    scc_lyot_stop,
    transmitted_lyot_field,
    reflected_lyot_field,
    psf_contrast,
    simulate_cred2,
    cdi_otf

include("configuration.jl")
include("masks.jl")
include("prescription.jl")
include("prepared.jl")
include("llowfs.jl")
include("diagnostics.jl")
include("movie.jl")

end # module ProperOptics
