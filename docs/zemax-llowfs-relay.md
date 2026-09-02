# Zemax LLOWFS relay: recovered prescription and numerical use

## Purpose and scope

This document records what was recovered from the SPIDERS Zemax archive, how
those data were converted into the scalar LLOWFS model now owned by
`SPIDERSSim.ProperOptics`, and
what remains uncertain. It covers the configuration-5 `LLOWFS main` path from
the reflective Lyot plane to the GoldEye detector. It also records the current
grid-convergence and CPU timing evidence because those results constrain how
the model should be used in an RTC.

The checked-in implementation is a centered, monochromatic, scalar Fresnel
model. It is not a full replacement for the Zemax sequential model. In
particular, it does not model polarization, coating phase or efficiency,
off-axis aberrations, surface figure, mechanical coordinates, or the physical
detector response.

## Data and provenance

The source is the user-supplied archive:

```text
SubaruPathFinderV11-oae-dm-beam-comp-f64-apo12-planar-customLens-coldStop-rightlocation-filter-narcisist-aperture.ZAR
SHA-256: 29ef0fe4fdb7a60931648700115b5fce3690604dadb4f9658a680b8957977666
```

The archive remains untracked because it may contain private instrument data.
The maintained local evidence copy is:

```text
/mnt/datadrive/DATA/SPIDERS/provenance/SubaruPathFinderV11-oae-dm-beam-comp-f64-apo12-planar-customLens-coldStop-rightlocation-filter-narcisist-aperture.ZAR
```

The prescription, catalog coefficients, derived values, and their provenance
are documented here so the implementation can be reviewed without publishing
the archive.

`zmxtools` was used to unpack the archive. The recovered `.zmx` file and the
embedded OHARA and fused-silica glass catalogs were inspected directly.
`ray-optics` was useful as an independent paraxial check, but its import did
not preserve the Zemax multi-configuration rules or square apertures. It was
therefore not treated as the prescription authority.

The Zemax file defines eight configurations. Configuration 5 is named
`LLOWFS main`; configuration 6 is `LLOWFS ref`. Configuration 5 ignores the
ordinary transmitted-Lyot surfaces 56--58 and follows the alternate LLOWFS
branch at surfaces 130--139. The Julia model uses configuration 5.

## Recovered optical path

| Zemax surface | Recovered element | Value used by `SPIDERSSim.ProperOptics` |
| --- | --- | --- |
| 130 | Lyot mirror | ideal flat mirror; reflective field comes from the complementary SCC Lyot mask |
| 131 | coordinate break and air | -5 degree fold, then 515 mm |
| 132 | doublet front | curvature +0.0076753039420361041 mm^-1, S-LAH71, 5.22416 mm |
| 133 | cemented interface | curvature -0.022355360815523563 mm^-1, S-NPH5, 1.8 mm |
| 134 | flat doublet back and air | 150 mm |
| 135 | filter | 2 mm F_SILICA |
| 136--137 | air to detector plane | 71.8759395036 mm; 3.84 x 4.8 mm square aperture |

Surface 130 contains `GLAS MIRROR`, but no patterned mask geometry or measured
coating prescription. The reflected LLOWFS field is therefore generated from
the complement of the modeled SCC Lyot-stop amplitude transmission. The
implementation uses amplitude coefficients whose squared magnitudes conserve
power at antialiased mask boundaries. A measured reflection phase, coating
efficiency, or figure map must replace this idealization when those data
become available.

The -5 degree coordinate break redirects the local optical axis. It is omitted
from the centered scalar propagation; this is valid for a model expressed in
the folded local coordinate system, but it does not reproduce global image
orientation or off-axis aberrations.

The Zemax Sellmeier and thermal coefficients are encoded in `src/llowfs.jl`.
At 1550 nm and 20 degrees C they give:

| Material | Refractive index |
| --- | ---: |
| S-LAH71 | 1.8134701603020393 |
| S-NPH5 | 1.8090114983950631 |
| F_SILICA | 1.4440236216697244 |

The corresponding first-order checks are:

| Quantity | Result |
| --- | ---: |
| Doublet effective focal length | 157.691076323 mm |
| Doublet back focal length | 153.859769413 mm |
| Lyot conjugate from the final doublet surface | 223.444722793 mm |
| Reduced physical detector distance | 223.260958398 mm |
| Detector displacement from the paraxial conjugate | -0.183764395 mm |
| Lateral magnification | -0.441273881834 |

The effective focal length varies from 157.513418 mm at 850 nm to
157.607199 mm at 1265 nm and 157.719631 mm at 1650 nm.
`SPIDERSSim.ProperOptics` computes
the indices at the wavelength of each call rather than fixing the 1550 nm
values.

Propagation through each glass thickness uses reduced distance `t / n` with
the vacuum wavelength. This is the paraxial scalar equivalent used by the
PROPER implementation. The curved air/glass and cemented interfaces retain
their paraxial optical powers; the flat rear and filter surfaces add only
reduced propagation distance.

Code locations:

- `src/proper_optics/llowfs.jl`: prescription, catalog dispersion, ABCD checks,
  and relay
  propagation;
- `src/proper_optics/prepared.jl`: prepared common-path and SCC propagation;
- `test/proper_optics.jl`: refractive-index, first-order, energy, sampling, and
  zero-allocation checks;
- the archived SpidersProper revision records the original standalone CPU and
  CUDA benchmark entry points.

## Detector and ROI evidence

The 4.8 x 3.84 mm Zemax aperture is exactly consistent with a 320 x 256 array
of 15 micrometre pixels. This is strong geometric evidence for the native
GoldEye format, but it is not a detector data sheet.

The deployed `spiders-working/imgadaptorservice` wire decoder reshapes camera
bytes as `(SizeX, SizeY)` and exposes `OffsetX` and `OffsetY`. It converts the
incoming integer pixels to `Float32` pixel-for-pixel. Inspection of the image
adaptor, calibration, and LOWFS paths found no binning, interpolation, or
resize stage. The measured calibration products are 34 x 47 pixels, so the
best current interpretation is a native-pixel hardware ROI or a software
crop. In VENOMS storage order, full frame is `(320, 256)` and the measured ROI
is `(34, 47)`. A conventional row-column optical array must map its axes
explicitly at this boundary.

Saved calibration FITS files do not retain the ROI origin or sufficient
orientation metadata. `OffsetX`, `OffsetY`, and an orientation/handedness
check must be captured from a live raw VENOMS frame before the simulated ROI
can be registered to the instrument. A centered ROI is only a provisional HIL
default.

Pixel integration, flux scaling, quantum efficiency, read noise, dark signal,
flat-field response, saturation, ADC behavior, exposure timing, and ROI
readout are intentionally outside this optical prescription. They belong in
the `AdaptiveOpticsSim` camera model and the HIL adapter.

## Grid selection and measured sensitivity

Matching one PROPER output sample to the native 15 micrometre detector pitch
at 1550 nm requires these beam-diameter fractions:

| Grid | Beam-diameter fraction |
| ---: | ---: |
| 256 | 0.4565295348866136 |
| 512 | 0.2282647674433068 |
| 1024 | 0.1141323837216534 |

The small detector ROI does **not** by itself permit the whole SPIDERS model
to use a small grid. Before the Lyot reflection, the numerical grid must still
sample the entrance pupil, focal-plane mask, SCC reference pinhole geometry,
Lyot-mask edges, and Fresnel padding. Cropping only the final intensity cannot
recover information lost or aliased in that common path.

Exploratory comparisons held the 15 micrometre detector sampling and pupil
sample count approximately fixed, then compared the centered normalized
47 x 34 ROI:

| Comparison | Normalized correlation | Relative L2 error |
| --- | ---: | ---: |
| 256 versus 1024 | about 0.785 | not retained |
| 512 versus 1024 | 0.97572 | 0.31085 |
| 1024 versus 2048 | 0.98745 | 0.0923 |

These results were obtained in exploratory runs; the raw comparison arrays
were not retained as a checked-in artifact. They show that 256 is not a
credible high-fidelity grid and that 512 is visually/structurally correlated
but not quantitatively converged. A 1024 grid is the current practical
reference, not a demonstrated asymptotic solution, and should still be
checked against 2048 for each intended operating regime.

## CPU timing evidence

The warmed prepared LLOWFS path was measured on an AMD Ryzen 7 6800H with
Julia 1.12.7, one Julia thread, one FFTW thread, `Float32`, and a 128 x 128 AO
residual. Preparation and compilation were outside the timed interval.

| Grid | Median | P95 | Mean frame rate | Warmed heap allocation |
| ---: | ---: | ---: | ---: | ---: |
| 256 | 37.780 ms | 38.906 ms | 26.38 frames/s | 0 bytes/call |
| 512 | 281.604 ms | 285.266 ms | 3.543 frames/s | 0 bytes/call |
| 1024 | 1182.104 ms | not measured with enough samples | about 0.846 frames/s | 0 bytes/call |

The comparable prepared SCC path at 256 had a 33.845 ms median. Thus most of
the 256-grid LLOWFS time is already spent in the common SPIDERS propagation;
making only the post-Lyot relay smaller cannot turn this model into a 1 kHz
plant.

Reproduce the CPU benchmark with, for example:

```bash
julia --project=. benchmarks/benchmark_prepared_llowfs.jl 512 100 5 128 Float32 0.2282647674433068
```

The benchmark prints environment metadata, raw latency samples, output
sampling, checksum, and warmed allocations. Raw timing output should be saved
with future performance claims.

## CUDA timing evidence

The CUDA path was measured on 2026-08-31 at source revision `48127e6` using an
NVIDIA GeForce RTX 3050 Ti Laptop GPU with 4 GiB, compute capability 8.6,
driver 596.08, CUDA.jl 6.3.1, CUDA runtime 13.3, Julia 1.12.6, one Julia
thread, `Float32`, and a 128 x 128 AO residual. The model used registered
Proper 0.2.0. CPU preparation, device adaptation, compilation, and ten warmup
calls were outside the samples.

The measured boundary is one closed-loop, synchronized call from the input
already resident on the GPU through the reflected-Lyot relay to sensor-plane
intensity. It includes kernel launches and `CUDA.synchronize()`. It excludes
host/device input transfer, detector simulation, AdaptiveOpticsSim Camera and
DM work, transport, queues, and RTC scheduling.

| Grid | Samples | Median | P95 | Maximum | Mean frame rate | Relative L2 versus CPU |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 200 | 3.316 ms | 4.138 ms | 5.105 ms | 293.813 frames/s | 3.04e-6 |
| 1024 | 100 | 13.773 ms | 14.195 ms | 14.922 ms | 72.261 frames/s | 3.30e-6 |
| 2048 | 30 | 54.667 ms | 56.443 ms | 56.478 ms | 18.261 frames/s | 3.27e-6 |

The 2048 P95 is only indicative because it has 30 samples. At the common 512
and 1024 grids, the CUDA median was about 85 times faster than the recorded
Ryzen 7 6800H CPU median. This is a hardware-specific comparison, not a general
CPU/GPU speedup claim.

An additional 256-grid run was both scientifically unsuitable according to
the convergence evidence and temporally noisy (7.247 ms median, 16.311 ms
P95, 114.145 mean frames/s). It is not used as evidence of CUDA scaling or RTC
capacity. The laptop GPU is shared with the Windows display and changed power
states automatically; no affinity, power locking, or exclusive-compute mode
was applied.

The synchronized CUDA call allocated about 111.8 KiB of host heap per 1024 or
2048 residual-input propagation. The HIL-style external-OPD input allocated
129,440 bytes at 1024. Allocation isolation attributed only 16 bytes to
`CUDA.synchronize()`; the remainder was inside `propagate_llowfs!`. The CPU
prepared path remains zero-allocation.

An allocation profile with `sample_rate=1.0` recorded 2,207 host allocations
and 117,520 sampled bytes for one warmed external-OPD call. The authoritative
`@allocated` total for that same boundary was 129,424 bytes before
synchronization. The largest recurring groups were CUDA/KernelAbstractions
submission objects created at Proper kernel launches:

- 24 quadratic-phase launches constructed `KernelCall`, `CompilerJob`,
  compiler-metadata, tuple, and managed-argument objects;
- 21 field-scaling launches constructed the same classes of host objects;
- GPU broadcast, centered `shift_copy!`, cubic-convolution resampling, and
  phase/mask operations contributed additional launch descriptors and views.

These are host submission allocations, not new detector outputs or large
device arrays. Removing them one Proper kernel wrapper at a time would require
changes across CUDA.jl/KernelAbstractions submission paths. CUDA Graph capture
provides a more direct amortization boundary for this fixed prepared model.

Reproduce the 1024 measurement with:

```bash
julia --project=. benchmarks/benchmark_cuda_prepared_llowfs.jl \
    1024 100 10 128 0.1141323837216534
```

These results support a roughly 72 Hz correctness-first 1024 optical plant on
this GPU, or about 18 Hz at 2048, before HIL integration overhead. They do not
support a greater-than-1-kHz full-PROPER LLOWFS plant.

### CUDA Graph experiment

The complete warmed propagation, including its cuFFT operations, was captured
once with Proper 0.3's general `prepare_cuda_graph` facility and then launched
repeatedly with `launch_cuda_graph!`. Capture and instantiation remained outside
the timed loop. This
reduced synchronized host allocation to 16 bytes per call; `CUDA.launch`
itself added zero observed host bytes, and the remaining 16 bytes came from
the benchmark's explicit synchronization.

The production-like experiments used `ExternalSpidersEntrance`, the same
prepared-input form intended for HIL. After capture, the existing device OPD
buffer was changed in place by 10 nm. Both graphs observed the new contents.
The LLOWFS and SCC outputs agreed with newly evaluated CPU references to
2.96e-6 and 5.14e-6 relative L2, respectively, confirming that capture did not
freeze the input values.

| Prescription/grid/input | Direct median | Graph median | Direct P95 | Graph P95 | Direct host bytes | Graph host bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LLOWFS, 512 residual | 2.994--3.096 ms | 2.752 ms | 3.736--4.048 ms | 2.872 ms | 111,440 | 16 |
| LLOWFS, 1024 external OPD | 14.307--14.368 ms | 14.141 ms | 14.544--14.926 ms | 14.409 ms | 129,440 | 16 |
| SCC, 1024 external OPD | 13.713--13.737 ms | 13.587 ms | 13.915--14.056 ms | 13.749 ms | 125,648 | 16 |

Each direct range brackets the graph phase with 100 samples per phase at 1024
and 200 per phase at 512. At 1024, device FFT work dominates, so graph capture
improved median latency by about one percent for both prescriptions while
eliminating almost all host allocation. Their measured graph throughputs were
70.6 LLOWFS frames/s and 73.6 SCC frames/s. At 512, launch overhead was more
important: graph capture improved LLOWFS median latency by roughly 8--11
percent and mean throughput by roughly 13--18 percent, with visibly tighter
tails in that run.

The graph is reusable only while the prepared model's array addresses, grid,
FFT plans, control flow, and stream-compatible operations remain stable.
Input *contents* may change in place. A new prepared object, grid size,
wavelength-dependent prescription, buffer replacement, or incompatible
configuration change requires recapture. Host-to-device transfer was not
inside the experiment; a CPU-generated DM/OPD frame must still be copied into
the stable device buffer, while a GPU-resident AdaptiveOpticsSim DM result can
feed it directly.

Standalone SPIDERS and the HIL loop use the same Proper graph executor. It can
capture the complete reflective Lyot `propagate_llowfs!` path or the complete
transmissive Lyot `propagate_scc!` path through the SCC focus. CUDA.jl's
high-level `@captured` macro captures and updates a graph on each macro
evaluation, which is convenient for changing operations but is not the
desired fixed hot-path lifecycle here.

The original standalone benchmark commands, retained in the archived source
revision, were:

```bash
julia --project=. benchmarks/profile_cuda_prepared_llowfs.jl \
    1024 0.1141323837216534 external
julia --project=. benchmarks/benchmark_cuda_graph_llowfs.jl \
    1024 100 0.1141323837216534 external
julia --project=. benchmarks/benchmark_cuda_graph_scc.jl \
    1024 100 0.1141323837216534 external
```

## Recommended two-grid LLOWFS method

Resampling is worth testing, but only at a physically and numerically safe
boundary: the complex reflected field at the Lyot plane. The proposed method
is:

1. Propagate the shared coronagraph on a grid established by convergence of
   the focal-plane and Lyot-plane fields.
2. Apply the reflective Lyot amplitude coefficient on that grid.
3. Crop a window containing all reflected complex-field support plus a guard
   band large enough to suppress wraparound.
4. Regrid the **complex field** (amplitude and phase) onto an LLOWFS-specific
   grid while preserving physical sampling, field center, and integrated
   power.
5. Propagate the recovered doublet/filter relay on the smaller grid.
6. Integrate onto native 15 micrometre pixels and then select the hardware or
   software ROI. Do not resize the detector intensity to make it fit.

This decouples the relay grid from the upstream coronagraph grid and can save
the FFT work unique to the LLOWFS branch. It does not reduce the cost of the
shared path, which currently dominates the CPU timing. A scaled Fresnel,
chirp-Z, or matrix-Fourier-style ROI propagator may eventually avoid computing
unused detector samples, but it must reproduce the slightly defocused relay;
a focal-plane-only transform is not sufficient.

The high-fidelity reference should retain the single large grid. The two-grid
path should be an explicitly selected prepared approximation, never an
implicit change to the reference result. For a greater-than-1-kHz RTC, the
likely online plant is a reduced LLOWFS response/Jacobian or another surrogate
generated and periodically checked by the high-fidelity PROPER model. Even a
successful post-Lyot regrid is unlikely to make the complete Fresnel plant run
at the requested rate.

### Correctness-first operating decision

On 2026-08-31, the user selected correctness over frame rate for the initial
integration. The initial HIL path will therefore use the existing single-grid
prepared propagation without a post-Lyot regrid. A 1024 grid is the current
practical working reference; selected 2048 runs remain necessary to measure
its error for the wavefront and DM-probe regimes used by the RTC. No real-time
frame-rate claim follows from this choice.

The two-grid method remains a later optimization. It will be introduced only
as an explicit alternative after the qualification below, with the
single-grid result retained as its comparison. This sequence avoids making
interpolation error part of the initial optical-model validation.

## Qualification required before enabling resampling

A candidate regrid must be compared with the single-grid 1024 and selected
2048 references over more than a flat wavefront. The comparison set should
include positive and negative tip, tilt, focus, and astigmatism probes; the
intended DM influence functions; and representative AO residuals. Record:

- complex-field energy before and after regridding;
- detector ROI correlation and relative L2 error;
- centroid, orientation, and pupil-image scale;
- finite-difference LLOWFS response/Jacobian error, including sign;
- sensitivity to crop width, guard band, relay-grid size, interpolation
  kernel, wavelength, and detector pixel integration;
- CPU and CUDA latency for the common path, regrid, relay, and detector
  sampling as separate stages.

Image similarity alone is not a sufficient RTC acceptance test: low-order
response error can matter even when two normalized pupil images look alike.
Numerical acceptance limits have not yet been agreed, so none are asserted
here.

## Confirmed assumptions

The user confirmed all five statements on 2026-08-31. They remain listed
because future instrument evidence may supersede them, and rejecting any one
would change the model or its deployment:

1. Configuration 5 (`LLOWFS main`), rather than configuration 6
   (`LLOWFS ref`), is the physical science path to simulate.
2. The 34 x 47 live image is a native-pixel ROI/crop with no binning or
   resampling between camera and LOWFS processing.
3. The online greater-than-1-kHz LLOWFS plant may be a reduced model, while
   full PROPER runs asynchronously for truth-model and HIL validation.
4. The two-grid approximation will be accepted using low-order response or
   Jacobian error, not normalized image correlation alone.
5. A live VENOMS header and an asymmetric calibration probe can be captured
   to establish ROI offsets, axis order, handedness, and sign.

## Gaps and unresolved questions

- The reflective Lyot mask has no recovered measured coating, phase, or
  surface-figure map.
- Zemax has not been run directly to export a complex field or detector image
  for end-to-end comparison with the Julia prescription.
- The ROI origin and detector orientation are unknown from saved products.
- The detector spectral response, gain, exposure behavior, noise, saturation,
  and pixel-response nonuniformity have not been identified here.
- The scalar model omits the -5 degree fold in global coordinates and cannot
  predict orientation or off-axis aberrations from that fold.
- CUDA timing has only been measured on one shared laptop GPU and has not been
  tested under sustained HIL load, fixed GPU power policy, or independent
  repetitions.
- CUDA Graph capture has been integrated through Proper 0.3 and
  AdaptiveOpticsProperHIL, and a changed in-place OPD buffer was verified. A
  sustained HIL loop with changing live DM frames has not yet been exercised.
- No two-grid regrid has yet been implemented or qualified.

The Zemax-to-Julia implementation and its unit checks first appeared in Git
commit `2c6b80c139ded488faecea59f5fe01eb33d30a11` (`Model Zemax LLOWFS relay`).
