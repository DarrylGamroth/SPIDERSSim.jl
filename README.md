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

