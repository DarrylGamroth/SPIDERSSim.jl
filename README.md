# SpidersProper.jl

`SpidersProper` is a Julia/Proper.jl port of the MATLAB end-to-end Fresnel
model for SPIDERS, the Subaru Pathfinder Instrument for Detecting Exoplanets
and Retrieving Spectra.

The port models the monochromatic Self-Coherent Camera (SCC) path:

```text
Subaru pupil -> f/13.901 input -> OAE1 -> DM -> OAE2 -> apodizer
             -> Tilt-Gaussian FPM -> Lens 0 -> SCC Lyot stop
             -> chopper plane -> Lens 1 -> SCC focal plane
```

It uses the pupil, apodizer, AO-residual, and focal-plane-mask FITS assets that
were supplied with the MATLAB model. The original `.m` files remain in the
repository as provenance and comparison references.

## Install

Proper `v0.2.0` is available from `DarrylGamrothRegistry`. Until SpidersProper
itself is registered, add the registry and install this repository explicitly:

```julia
using Pkg
Pkg.Registry.add(Pkg.RegistrySpec(
    url="https://github.com/DarrylGamroth/PackageRegistry.git",
))
Pkg.add(url="https://github.com/DarrylGamroth/SpidersProper.jl.git")
```

For local development, use `Pkg.develop(path=...)` instead.

## Run the prescription

```julia
using SpidersProper

config = SpidersConfig(reference_pinhole=true)
result = spiders_proper(1.25e-6, 1024; config)

intensity = spiders_intensity(result)
field = spiders_field(result)

result.wavefront.sampling_m
result.lyot_stop.main_pupil_transmission
result.lyot_stop.pinhole_transmission
result.aperture_throughput
```

The `reference_pinhole` default is `false`, matching `phFlag = 0` in
`spidersProper5.m`. Enable it for an SCC reference beam, as above.

The bundled FITS pupil and apodizer are used by default. An analytic Subaru
pupil and the original radial-polynomial apodizer are also available:

```julia
config = SpidersConfig(
    pupil_mode=:analytic,
    apodizer_mode=:polynomial,
    fpm_type=:TGV,
    reference_pinhole=true,
)
result = spiders_proper(1.55e-6, 1024; config)
```

For the supplied 1,250 nm FPM phase map:

```julia
config = SpidersConfig(
    fpm_mode=:map,
    fpm_map_path=joinpath(pkgdir(SpidersProper), "FPM-J-1250nm-Hlevels-500nmPixel.fits"),
    fpm_map_sampling_m=0.5e-6,
    fpm_map_type=:phase,
)
```

Measured mirror maps can be enabled by providing their paths:

```julia
config = SpidersConfig(
    oae1_error_path="/path/to/oae1v2.6SagMicronProfilo.fits",
    oae2_error_path="/path/to/oae2v2FigureErrorMicronProfiloTTFAremoved.fits",
    dm_correction=true,
)
```

## AO movie cubes

`spiders_movie` accepts AO residuals directly instead of depending on the
original hard-coded PASSATA `.mat` path:

```julia
# residuals_nm has dimensions (y, x, time) and values in nanometers.
cube = spiders_movie(
    1.55e-6,
    residuals_nm;
    gridsize=1024,
    config=SpidersConfig(reference_pinhole=true),
)
```

The output is allocated as a dense `Float64` cube. A 1024 x 1024 x 1000 run
therefore requires about 7.8 GiB for the output alone.

## Diagnostics

- `psf_contrast` returns radial 25th/50th/75th-percentile profiles in
  `lambda/D`. Pass the peak of a non-coronagraphic PSF as `reference_peak` to
  obtain normalized contrast.
- `simulate_cred2` implements the original simplified photon/read-noise camera
  model without plotting.
- `cdi_otf` returns the SCC OTF decomposition, visibility map, and the two CDI
  planet-intensity estimates.

## Fidelity and deliberate changes

- The physical spacings, focal lengths, clear apertures, mask sizes, and
  propagation order follow `spidersProper5.m`.
- Plotting was separated from propagation; the Julia prescription has no
  unconditional figure creation or global colormap changes.
- The self-contained defaults disable OAE surface errors and DM correction,
  because the referenced metrology FITS files were not supplied.
- The bundled pupil/apodizer FITS maps replace the MATLAB default's missing
  `prop_subaruPupilSpiders` helper. Select `pupil_mode=:analytic` to use the
  ported geometric approximation.
- The model remains monochromatic per call. It does not add the LLOWFS,
  calibration-source, iFTS, coating, polarization, or optomechanical paths.
- The detector helper retains the MATLAB model's simplified throughput and
  zero-point assumptions. For large Poisson means it uses the Gaussian limit.

## Test

```julia
using Pkg
Pkg.test()
```

The tests exercise both pupil/apodizer implementations, TG/TGV masks, Lyot
stop construction, an end-to-end SCC propagation, detector/CDI diagnostics,
and an in-memory AO movie frame.

## Propagation benchmark

The standalone benchmark separates the first propagation after package import
from warmed, closed-loop propagation samples. Its operation boundary is the
current public `spiders_proper` call, including repeated FITS reads, mask and
workspace construction, and the returned wavefront:

```bash
julia --project=. benchmarks/benchmark_propagation.jl 512 20 3 128
```

The arguments are propagation grid size, sample count, warmup count, and AO
residual grid size. Use at least 100 samples before interpreting the reported
p99. The script retains the raw latency samples in its output and performs a
separate allocation probe.

For the prepared repeated-frame path:

```bash
julia --project=. benchmarks/benchmark_prepared_propagation.jl 512 20 3 128
```

`prepare_spiders` loads static FITS data, builds focal/Lyot and clear-aperture
masks, and allocates PROPER workspaces once. `spiders_propagate!` then reuses
the bound input, wavefront, scratch, masks, FFT plans, and output intensity:

```julia
residual_nm = zeros(Float32, 128, 128)
prepared = prepare_spiders(
    1.25e-6,
    512;
    config=SpidersConfig(reference_pinhole=true),
    ao_residual_nm=residual_nm,
)

residual_nm .= next_residual
intensity = spiders_propagate!(prepared)
```

The warmed CPU operation has a zero-byte Julia heap-allocation contract on the
tested Julia/Proper versions. Inputs and the returned output are bound mutable
buffers: update their contents between frames and do not replace them.

## AdaptiveOpticsSim integration boundary

AdaptiveOpticsSim should own the atmosphere, telescope pupil, RTC commands,
and physical DM response. Bind its caller-owned science-pupil amplitude and
OPD storage directly to the prepared SPIDERS path:

```julia
prepared = prepare_spiders(
    wavelength_m,
    propagation_size;
    config,
    pupil_amplitude=pupil_amplitude(science_pupil),
    pupil_sampling_m=telescope_diameter_m / pupil_resolution,
    opd_m=opd_map(science_pupil),
    T=Float32,
)

# Each model-time step updates the already-bound pupil arrays.
update_surface!(dm)
apply_surface!(science_pupil, dm, DMAdditive())
relative_intensity = spiders_propagate!(prepared)
```

This is the collapsed residual-OPD seam: the atmosphere, telescope errors,
and DM OPD are composed on the entrance-pupil grid before PROPER propagation.
Do not also enable `SpidersConfig.dm_correction`; the prepared path rejects it
to prevent applying the DM twice. If the as-built SPIDERS DM conjugate must be
modeled explicitly, AdaptiveOpticsSim should still own commands and influence
functions, but a future prepared input should apply its resulting OPD at the
prescription's DM plane rather than collapsing it into the entrance pupil.

The final SPIDERS array is dimensionless, cell-integrated relative intensity.
Before passing it to an AdaptiveOpticsSim detector, prepare an explicit,
flux-conserving mapping from the PROPER focal sampling to the physical detector
pixel pitch and declare an entrance-photon-rate scale. Source magnitude,
collecting area, bandpass, and optical throughput belong in that radiometric
scale, not in the detector.

For a C-RED2-like camera, prefer AdaptiveOpticsSim's `Detector` with an
`InGaAsSensor`. The local `simulate_cred2` helper is retained for MATLAB parity
and combines several concerns: focal-grid resampling, a fixed stellar zero
point and throughput, Poisson noise, Gaussian read noise, gain conversion, and
integer flooring. AdaptiveOpticsSim separates those concerns and additionally
supports calibrated QE, dark current, glow, persistence, nonlinearity,
saturation, response/defect maps, readout timing and modes, quantization, and
prepared zero-allocation acquisition. A parity starting point is photon plus
read noise at 25 electrons RMS, QE 0.7, zero dark current, and an InGaAs sensor;
use measured camera values before treating it as a C-RED2 qualification.

Use `AdaptiveOpticsProperHIL.jl` for the integration-level benchmark. That
benchmark should time the prepared AO-command-to-SPIDERS-observation boundary,
with caller-owned OPD and pupil buffers. The standalone result is a baseline
for the unprepared prescription and is not a HIL latency claim.

## Package and deployment roles

Keep this optical prescription as an independent `SpidersProper.jl` package.
It owns the SPIDERS optical layout, calibration assets, prepared propagation,
and MATLAB-parity helpers, and depends on `Proper.jl`.

Keep `AdaptiveOpticsSim.jl` independent of PROPER and instrument-specific
prescriptions. It owns the atmosphere, telescope, wavefront-sensor and detector
physics, DM commands and influence functions, and RTC simulation.

Use `AdaptiveOpticsProperHIL.jl` as the deployable integration environment. It
should depend on `AdaptiveOpticsSim`, `Proper`, and `SpidersProper`, bind the
AdaptiveOpticsSim pupil/OPD arrays into `PreparedSpiders`, map the propagated
focal plane into an AdaptiveOpticsSim detector, and own HIL profiles,
benchmarks, process wiring, and cadence. Commit that application's Manifest so
a deployed HIL configuration pins exact dependency revisions; do not commit a
Manifest in this library package.

Recommended release order:

1. Use the registered `Proper v0.2.0`, which contains the prepared API.
2. Tag `SpidersProper v0.1.0` and add it to `DarrylGamrothRegistry`.
3. Add `SpidersProper` to `AdaptiveOpticsProperHIL.jl`, commit its Manifest,
   and validate the complete prepared propagation and detector path there.
