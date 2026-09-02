# SPIDERS Simulation Agent Instructions

This package owns the SPIDERS instrument simulation. It depends on
AdaptiveOpticsSim.jl for reusable AO physics and graph execution,
AdaptiveOpticsProperHIL.jl for the reusable AOS–Proper graph adapter, and
Proper.jl for physical-optics propagation.

- Keep SPIDERS configuration, optical prescriptions, detector mappings,
  chopper behavior, calibration provenance, graph files, RTC integration, and
  instrument validation here.
- Do not duplicate reusable AOS algorithms or the generic Proper graph adapter.
- Keep inferred and placeholder values explicitly distinguished from measured
  or deployment-configured values.
- Separate immutable optical/calibration plans, persistent instrument state,
  replaceable workspace scratch, and caller-visible products.
- Preserve the deployed `(SizeX, SizeY)` image convention and exact actuator
  ordering at external boundaries.
- Use deterministic inputs and explicit seeded RNG ownership in validation.
- Keep private `/mnt/datadrive/DATA` products as execution evidence; committed
  tests must use synthetic or redistributable fixtures.
- Validate CPU first. GPU validation is optional and belongs in explicit test
  or benchmark environments.
- Do not commit root manifests unless this package is deliberately converted
  from a reusable package into an application environment.

