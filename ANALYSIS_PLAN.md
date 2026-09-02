# SPIDERS Simulation Analysis Plan

## Objective

Create a reusable SPIDERS instrument simulation that produces physically and
operationally traceable LLOWFS and SCC detector products, supports the
camera-synchronous SCC chopper, and can validate SPIDERS control algorithms in
an independent RTC process. Success means the package can distinguish optical
model defects, calibration defects, and RTC integration defects with
deterministic evidence.

## Data and provenance

- `~/workspaces/spidersProperPackage`: Julia/MATLAB Proper prescriptions,
  Zemax-derived relay notes, masks, and optical-development fixtures.
- `~/workspaces/codex/spiders` (also `RTC_planning/spiders`): deployed legacy
  LLOWFS, SCC, chopped-frame sorting, calibration, integrator, and DM-mixer
  behavior. These repositories are requirements/parity references, not code to
  copy without review.
- `SpiderMan.jl` and the deployed SPIDERS service configuration: mechanical
  state, stage presets, wheel serials, and slot tables.
- `/mnt/datadrive/DATA/LOWFS` and `/mnt/datadrive/DATA/SCC`: private calibration
  arrays and observations. Reconnaissance records paths, shapes, element types,
  dates, and checksums without committing private arrays.
- Proc. SPIE 12185, DOI `10.1117/12.2629644`: documented regular-H optical
  design. Paper values do not supersede as-built calibration.

Known limitations: the installed reference pinhole, wheel-to-arm ownership,
qualified detector sampling/crop/orientation, as-built focal-plane-mask map,
apodizer map, and several selector states remain unconfirmed. Camera ownership
is known: LLOWFS uses GoldEye and SCC uses C-RED 2.

## Environment and maturity

- Language and environment: Julia 1.12 package; Proper.jl physical optics;
  AdaptiveOpticsSim algorithm graphs; optional Python pyRTC integration in an
  isolated test environment.
- Maturity: reusable package with explicitly provisional optical profiles.
- Existing packages: AdaptiveOpticsSim.jl, AdaptiveOpticsProperHIL.jl,
  Proper.jl, and SpidersProper.jl.

## Working stance

- SPIDERSSim owns the instrument; AdaptiveOpticsProperHIL owns only reusable
  cross-package integration.
- Correctness and traceability precede performance optimization.
- CPU is the qualification reference. GPU execution must preserve the same
  optical model and explicitly copy products at a host RTC boundary.
- LLOWFS and SCC have separate cadences and graph boundaries.
- A pyRTC implementation may reuse pyRTC lifecycle, shared-memory, and loop
  infrastructure, but its scientific processing must match reviewed SPIDERS
  behavior rather than masquerade as a Shack–Hartmann or Pyramid sensor.

## Target outputs

- Deterministic 34×47 LLOWFS and 416×380 SCC relative-intensity products.
- Deterministic fringed/unfringed SCC pairs and their phase-tagged difference.
- A qualification report for every provisional optical and mechanical choice.
- LLOWFS calibration products mapping normalized masked image residuals and
  centroid tip/tilt into the deployed modal/actuator basis.
- SCC calibration products mapping chopped differential images into the
  selected dark-hole actuator basis.
- Independent-process RTC fixtures that accept complete detector products and
  return exact-sequence actuator commands.
- Human-inspection viewers for raw, reference, residual, chopped, and corrected
  products.

## Working knowledge

- 2026-09-01: deployed products use `(SizeX, SizeY)` array order. The maintained
  provisional shapes are LLOWFS `(34, 47)` and SCC `(416, 380)`.
- 2026-09-01: LLOWFS is the GoldEye arm and SCC is the C-RED 2 arm. Deployed
  dark/reference artifacts independently confirm Float32 products with shapes
  `(34, 47)` and `(416, 380)`, respectively. C-RED 2 calibration scripts also
  establish a 640×512 full detector coordinate system.
- 2026-09-01: the legacy LLOWFS path normalizes selected pixels by total flux,
  subtracts a reference, projects image residuals into modes, independently
  estimates centroid tip/tilt, and maps modes to actuator coefficients.
- 2026-09-01: the legacy SCC path subtracts a reference on a selected dark-hole
  mask and applies an image-to-actuator matrix. The old implementation's
  phase/amplitude masking requires review before reuse.
- 2026-09-01: the legacy chopped-frame sorter infers phase from frame parity and
  running flux. The simulation already knows the commanded chopper phase and
  should publish it explicitly; parity/flux inference is a deployment-compatibility
  adapter, not the canonical model.

## Chunks

### SP-001: Establish package ownership
- Description: Move SPIDERS profiles, prescriptions, chopper state, graph
  files, and their CPU tests into SPIDERSSim.jl while retaining the generic
  Proper graph adapter in AdaptiveOpticsProperHIL.jl.
- Depends on: none
- Verification: both packages load and pass their focused CPU tests; no
  SPIDERS public symbol remains exported by AdaptiveOpticsProperHIL.
- Status: complete
- Notes: SPIDERSSim baseline `af3e5d5`; generic-adapter extraction
  `AdaptiveOpticsProperHIL.jl@c812ae5`.

### SP-002: Reconcile optical artifacts
- Description: Compare the provisional LLOWFS/SCC outputs with SpidersProper,
  Zemax notes, SpiderMan settings, and private calibration-product metadata.
- Depends on: SP-001
- Verification: deterministic inspection script plus recorded shapes,
  checksums, assumptions, and residual images.
- Status: not-started
- Notes: do not commit private observatory data.

### SP-003: Define detector and chopper products
- Description: Model GoldEye/LLOWFS and C-RED 2/SCC acquisition boundaries,
  exact ROI order, phase tags, dropped-frame behavior, and chopped difference
  ownership.
- Depends on: SP-001, SP-002
- Verification: synthetic phase sequence, reset, dropped-frame, orientation,
  and detector-response tests.
- Status: in-progress
- Notes: GoldEye and C-RED 2 ownership, Float32 product shapes, and deployed
  `(SizeX, SizeY)` order are now explicit. Phase-tagged rolling difference,
  sequence monotonicity, reset, and dropped-frame invalidation are implemented.
  The pyRTC viewer demo publishes both detector products and the adjacent
  chopped difference through the maintained shared-memory convention.
  Detector-response qualification remains.

### SP-004: Implement LLOWFS calibration and processor
- Description: Port reviewed normalization, reference subtraction, masked
  projection, centroid tip/tilt, modal-to-actuator mapping, and validity checks
  into typed prepared Julia and pyRTC-facing forms.
- Depends on: SP-002, SP-003
- Verification: analytical synthetic fixtures plus parity against selected
  deployed calibration artifacts.
- Status: in-progress
- Notes: the typed CPU processor now implements total-flux normalization,
  selected-pixel reference subtraction, modal reconstruction, centroid
  tip/tilt reporting, and actuator mapping. Artifact parity remains.

### SP-005: Implement SCC calibration and processor
- Description: Define chopped-difference formation, dark-hole masks,
  interaction-matrix generation, regularization, and actuator reconstruction.
- Depends on: SP-002, SP-003
- Verification: synthetic complex-speckle probes, push-pull linearity, retained
  rank, and parity against selected SCC artifacts.
- Status: in-progress
- Notes: the typed CPU processor now applies an explicit correction mask and
  maps a complete C-RED 2 chopped difference to actuator coefficients. Rolling
  pair formation is implemented; interaction-matrix generation and artifact
  parity remain.

### SP-006: Validate independent-process control
- Description: Add isolated pyRTC-compatible LLOWFS and SCC processes using
  POSIX shared memory, explicit sequence/phase metadata, and complete actuator
  commands.
- Depends on: SP-004, SP-005
- Verification: calibration followed by deterministic closed loops with
  residual-image, contrast, OPD, command, and failure-containment gates.
- Status: not-started
- Notes:

### SP-007: Qualify accelerator execution
- Description: Validate CPU/AMDGPU/CUDA optical equivalence and explicit GPU to
  host RTC handoff after the CPU model is accepted.
- Depends on: SP-002, SP-003
- Verification: numerical tolerances, warmed allocation evidence, and recorded
  frame-service benchmarks.
- Status: in-progress
- Notes: complete static LLOWFS and SCC optical graphs are capture-qualified.
  The 2026-09-01 baseline establishes CPU, AMDGPU/HIP Graph, and CUDA Graph
  numerical agreement and synchronized service times at the explicit host RTC
  boundary. CUDA capture reduces latency and host allocation on the measured
  RTX 3050 Ti; HIP capture is latency-neutral on the measured Radeon 680M but
  removes approximately 173--176 KiB of host allocation per graph step.
  Detector response, chopping, and closed-loop command application remain
  outside this benchmark and are not yet accelerator-qualified.

## Open questions

- Which AO3k/woofer/tweeter actuator surface should each LLOWFS and SCC command
  target, and what is the exact deployed command order and unit?
- Is the installed SCC reference pinhole 203 µm or 291 µm, and what is its
  pupil-plane position angle?
- Which science-wheel serial belongs to each arm, and which slots define the
  normal observing configuration?
- Should the first external RTC fixture run the legacy linear processors, the
  SPIDERS V2 service implementations, or both as separate parity targets?
- What contrast or residual metric is the authoritative SCC closed-loop gate?

## Evidence ledger

- 2026-09-01: existing AdaptiveOpticsProperHIL CPU tests establish finite,
  nonnegative, input-sensitive LLOWFS/SCC products and allocation-free warmed
  graph steps for the provisional prescriptions.
- 2026-09-01: existing chopper tests establish deterministic fringed/unfringed
  alternation, phase status aligned with graph sequence, reset semantics, and
  allocation-free warmed stepping.
- 2026-09-01: focused native processor tests establish deterministic GoldEye
  normalization/modal reconstruction and C-RED 2 masked difference
  reconstruction on dense and non-contiguous host views. Both warmed process
  calls allocate zero Julia heap bytes and pass with bounds checks enabled.
- 2026-09-01: C-RED 2 pairing tests establish explicit phase ownership,
  rolling fringed-minus-unfringed products, monotonic sequence rejection,
  dropped-frame invalidation, recovery on the next adjacent pair, reset, and
  zero-allocation warmed acceptance.
- 2026-09-01: clean-revision production-grid benchmarks establish relative L2
  agreement within `5.3e-6` between CPU and both accelerator backends. The
  raw 100-sample, three-repetition records and environment provenance are in
  `benchmark/results`.
- 2026-09-01: the AMDGPU SPIDERS pyRTC viewer completed 1,950 paced frames
  with exactly 975 fringed and 975 unfringed C-RED 2 exposures, live GoldEye
  and SCC products, a nonzero adjacent chopped difference, early viewer-close
  handling, and clean shared-memory teardown. This is visualization evidence,
  not an RTC closed-loop result.
