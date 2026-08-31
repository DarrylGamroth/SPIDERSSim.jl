module SpidersProper

using FFTW: fft, fftshift, ifft, ifftshift
using Proper
using Random: AbstractRNG, default_rng, rand, randn
using Statistics: median, quantile

include("config.jl")
include("masks.jl")
include("prescription.jl")
include("prepared.jl")
include("diagnostics.jl")
include("movie.jl")

export SpidersConfig, SpidersResult, LyotStop
export spiders_proper, spiders_proper5, spiders_movie
export PreparedSpiders, prepare_spiders, spiders_propagate!
export spiders_intensity, spiders_field
export subaru_pupil!, radial_apodizer, tilt_gaussian_fpm, scc_lyot_stop
export psf_contrast, simulate_cred2, cdi_otf

end
