# Optical graph baseline: 2026-09-01

This baseline measures one complete, synchronized static optical graph step at
closed-loop concurrency one. Both graphs use the provisional regular-H profile:
a 512 x 512 propagation grid, a 128 x 128 pupil, a 34 x 47 LLOWFS product, and
a 416 x 380 SCC product. Each phase contains 100 samples in each of three
independent repetitions after 10 warmup steps.

The `publish` phase adds a complete copy into preallocated host storage. It is
the current boundary for a CPU RTC. Compilation, FFT planning, graph capture,
detector noise, SCC chopping, RTC reconstruction, and DM updates are excluded.

| Device and mode | Arm | Core p50 (ms) | Core p99 (ms) | Mean fps | Publish p50 (ms) | Host bytes/step |
|---|---:|---:|---:|---:|---:|---:|
| Ryzen 7 6800H CPU, direct | LLOWFS | 116.542 | 119.862 | 8.56 | 116.220 | 0 |
| Ryzen 7 6800H CPU, direct | SCC | 128.611 | 132.132 | 7.79 | 129.416 | 0 |
| Radeon 680M, HIP direct-after | LLOWFS | 3.980 | 4.433 | 247.88 | 3.989 | 173,104 |
| Radeon 680M, HIP Graph | LLOWFS | 3.934 | 4.587 | 251.87 | 3.932 | 0 |
| Radeon 680M, HIP direct-after | SCC | 4.670 | 4.824 | 213.97 | 4.783 | 175,520 |
| Radeon 680M, HIP Graph | SCC | 4.642 | 5.176 | 212.70 | 4.737 | 0 |
| RTX 3050 Ti, CUDA direct-after | LLOWFS | 1.763 | 2.273 | 547.08 | 1.721 | 123,904 |
| RTX 3050 Ti, CUDA Graph | LLOWFS | 1.327 | 1.424 | 754.28 | 1.360 | 0 |
| RTX 3050 Ti, CUDA direct-after | SCC | 2.509 | 3.303 | 390.48 | 2.306 | 125,808 |
| RTX 3050 Ti, CUDA Graph | SCC | 1.716 | 1.892 | 578.18 | 1.955 | 0 |

`direct-after` is the direct execution phase repeated after capture measurement;
the raw records also retain `direct-before` to expose drift and warm-clock
effects. Frames per second is the reciprocal of the arithmetic mean service
time, not the reciprocal of p50.

CUDA Graph reduced core p50 by 25% for LLOWFS and 32% for SCC relative to the
same run's direct-after phase. HIP Graph was approximately latency-neutral on
this workload, while both command-graph paths reduced measured Julia host
allocation to zero. Accelerator results agreed with the CPU reference within
`5.3e-6` relative L2 for both a flat wavefront and a changed `2e-7` m tilt.

The systems are not vendor-matched. The HIP result used an integrated Radeon
680M (`gfx1030`) with AMDGPU.jl 2.8.0, HIP 7.2, and rocFFT 1.0.36. The CUDA
result used an RTX 3050 Ti Laptop GPU with CUDA.jl 6.3.1, CUDA runtime 13.3,
and NVIDIA Windows/WSL driver 596.08. CPU and HIP ran with Julia 1.12.7 on
Linux; CUDA ran with Julia 1.12.6 under WSL2.

All three records use these clean source revisions:

- SPIDERSSim.jl `de3ef5a6ad5d4a8ba7fb680a3d4ec13fe29ba46e`
- AdaptiveOpticsProperHIL.jl `f0865ebdfd2c02f7ac2d6dd992ef8050dd233894`
- AdaptiveOpticsSim.jl `5c17514afb95b35e4c370bb9ccab67801bf302a4`
- Proper.jl `409171c1c3f12e869ac3b5484a269240dd4c0c7e`

The JSON files retain every raw latency sample, correctness result, repetition
summary, backend versions, kernel, thread counts, and source cleanliness.
