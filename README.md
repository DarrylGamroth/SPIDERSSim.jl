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
