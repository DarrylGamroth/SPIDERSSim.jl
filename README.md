# SPIDERSSim.jl

Instrument-level SPIDERS simulation built on AdaptiveOpticsSim.jl, Proper.jl,
and the reusable Proper graph adapter in AdaptiveOpticsProperHIL.jl.

The package owns SPIDERS optical configurations, LLOWFS and SCC prescriptions,
camera-synchronous SCC chopping, detector mappings, calibration provenance,
and external-RTC validation. Its current regular-H model is explicitly
provisional: it is useful for topology, software integration, and calibration
development but is not yet an as-built optical model.

The maintained implementation is being extracted from
AdaptiveOpticsProperHIL.jl. See [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md) for the
validation sequence and unresolved instrument choices.

The first native control processors are also available. `prepare_llowfs`
seals the GoldEye reference, masks, modal calibration, and actuator mapping;
`process_llowfs!` consumes one complete host frame without allocating after
warmup. `prepare_scc` seals the C-RED 2 chopped-difference reference, dark-hole
masks, and image-to-actuator calibration; `process_scc!` consumes the complete
fringed-minus-unfringed product. The processor contracts intentionally require
one-based, host-scalar-indexable matrices; GPU optical products cross an
explicit device-to-host RTC boundary first.

`prepare_scc_pairing` creates the persistent C-RED 2 acquisition state.
`accept_scc_frame!` consumes an explicitly phase-tagged complete frame and
publishes a rolling fringed-minus-unfringed product only when the two retained
phases are adjacent in sequence. A dropped frame invalidates the product until
the next adjacent pair, so stale images cannot silently cross a gap.

## BAX307 measured deformable mirror

The instrument package validates the ALPAO BAX307 command map, actuator-health
map, sampled support, and measured influence matrix before preparing the
reusable AOS measured-DM graph node. FITS loading is optional; load `FITSIO` to
activate it and state the source influence unit explicitly:

```julia
using FITSIO
using SPIDERSSim

calibration = load_bax307_calibration(
    "/mnt/datadrive/DATA";
    influence_sample_unit_m=1e-6, # provisional until metrology confirms it
)
best_flat = load_bax307_best_flat("/mnt/datadrive/DATA")

dm = bax307_deformable_mirror_node(
    :bax307,
    calibration;
    pupil_diameter_m=0.033,
)
parameters = bax307_graph_parameters(:bax307, calibration) # host target
```

The plan snapshots host calibration arrays and preserves the exact 468-command
ordering. Its response model clips the complete normalized command and zeros
actuators marked invalid. Pass the selected AOS graph target as the third
argument to `bax307_graph_parameters` to prepare CUDA or AMDGPU startup
parameters. Calibration file paths and FITS semantics remain outside
AdaptiveOpticsSim core.

## Optical graph benchmarks

The maintained benchmark compares ordinary stream execution with complete
CUDA/HIP Graph capture for the production-sized provisional LLOWFS and SCC
optical graphs. It records synchronized core latency, complete device-to-host
publication latency, host allocations, raw samples, correctness against CPU,
source revisions, and the runtime environment:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/optical_graph.jl cpu 100 3
julia --project=benchmark benchmark/optical_graph.jl amdgpu 100 3
julia --project=benchmark benchmark/optical_graph.jl cuda 100 3
```

These are closed-loop service-time measurements with one outstanding frame;
they are not fixed-arrival-rate capacity claims. CUDA and HIP results are
device-specific and must not be compared as though they ran on identical
hardware. See [`benchmark/results`](benchmark/results/README.md) for the
recorded baseline and raw samples.

## Live pyRTC viewer

The SPIDERS viewer demo uses HR 8799 as an on-axis H-band point source. This is
a real [SCExAO/CHARIS observing target](https://subarutelescope.org/Projects/SCEXAO/scexaoWEB/080images.web/030science.web/indexm.html);
the adopted 2MASS H magnitude is
[5.280](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=HR+8799). A
deterministic, provisional four-layer atmosphere evolves above the Subaru
entrance pupil. The demo publishes masked pupil-plane atmospheric OPD, the
source-scaled GoldEye LLOWFS photon-arrival-rate frame, the most recent
unfringed C-RED 2 SCC photon-arrival-rate image, the paired
fringed-minus-unfringed SCC product, and the current measured BAX307
deformable-mirror surface OPD through pyRTC-compatible shared memory.

The atmospheric profile is a visualization fixture, not a reconstruction of a
specific HR 8799 observing night. The source scaling includes catalogued H-band
brightness and modeled geometric propagation, but not qualified instrument
throughput, exposure, or detector response. When `/mnt/datadrive/DATA` is
available, the viewer validates and loads the BAX307 24-by-24 actuator map,
465-actuator health map, 128-by-128 manufacturer influence product, and
468-element best flat. Set `SPIDERS_DATA_ROOT` to select another calibration
root or `SPIDERS_BAX307_COMMAND=zero` to start from zero command. The available
influence FITS product does not declare trustworthy physical units, so the
viewer currently makes the explicit provisional assumption
`SPIDERS_BAX307_INFLUENCE_SAMPLE_UNIT_M=1e-6` metres of mirror surface per
stored sample. Because the optical surrogate does not yet contain measured OAE
maps, the viewer balances the selected best-flat surface with an inferred
equal-and-opposite static OPD. This exposes the physical DM surface without
pretending that an otherwise absent static aberration is left uncorrected.
Neither that inferred aberration nor the influence units are qualified, so
they must be replaced before interpreting absolute contrast.

Set `PYRTC_PYTHON` to the Python interpreter in an environment containing the
GitHub pyRTC package and viewer:

```bash
julia --project=examples/pyrtc -e 'using Pkg; Pkg.instantiate()'
PYRTC_PYTHON=/path/to/venv/bin/python \
  julia --project=examples/pyrtc \
  examples/pyrtc/run_viewer_demo.jl amdgpu 120 20
```

Use `cpu` or `cuda` in place of `amdgpu` as appropriate. This demonstration
visualizes the provisional optical and camera-synchronous chopper products; it
does not yet run the LLOWFS or SCC reconstruction and control processes.
