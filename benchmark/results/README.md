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

## Shared two-arm propagation: 2026-09-02

The maintained `provisional_spiders_optical_node` propagates the common
entrance-pupil field through the focal-plane mask and into the Lyot plane once.
It then copies that coherent Lyot-plane field into the complementary LLOWFS
and SCC branches. Tests establish bitwise equality with two independent CPU
prescriptions for flat and tilted OPD inputs. The benchmark uses the chopped
SCC prescription and therefore measures the actual two-product optical
boundary used by the pyRTC viewer.

The pupil and DM grids remain 128 by 128 in every row. Only the padded Proper
propagation grid changes. Timings use 100 samples in each of three repetitions
after 10 warmup steps, one outstanding frame, one Julia thread, and one FFT
provider thread.

| Device | Propagation grid | Independent p50 (ms) | Shared p50 (ms) | Shared p99 (ms) | Shared mean fps | Shared publish p50 (ms) | Host bytes/step |
|---|---:|---:|---:|---:|---:|---:|---:|
| Ryzen 7 6800H CPU | 512 | 240.166 | 148.396 | 154.510 | 6.71 | 148.177 | 0 |
| Radeon 680M, HIP stream | 512 | 8.427 | 5.999 | 6.258 | 166.20 | 5.692 | 232,032 |
| Radeon 680M, HIP stream | 1024 | 45.078 | 28.219 | 29.202 | 35.35 | 29.352 | 234,544 |
| Radeon 680M, HIP stream | 2048 | 171.769 | 111.047 | 116.472 | 8.98 | 111.172 | 234,992 |
| RTX 3050 Ti, CUDA stream | 512 | 3.859 | 2.514 | 3.494 | 402.19 | 2.318 | 167,856 |
| RTX 3050 Ti, CUDA stream | 1024 | 13.903 | 9.019 | 10.698 | 109.14 | 9.234 | 170,112 |
| RTX 3050 Ti, CUDA stream | 2048 | 54.709 | 37.029 | 92.723 | 24.61 | 41.016 | 170,112 |

`publish` is a separately measured phase, not the core measurement plus an
arithmetically isolated copy time; normal run-to-run variation can therefore
make its median slightly lower. It copies both complete products to
preallocated host arrays. The shared node reduced core p50 by 29--38% over
duplicated propagation on these measured configurations.

This chopped workload currently uses stream execution. It is intentionally
not described as CUDA Graph or HIP Graph evidence: chopper phase adoption and
the shared wavefront fork have not yet been capture-qualified. The reported
GPU host allocations are therefore real remaining overhead, despite being
lower than the independent-node path. The earlier fixed-arm command-graph
baseline above measures a different boundary.

Every accelerator result passed a CPU reference gate for both detector
products. The largest retained relative L2 residual was `4.94e-5` at the 2048
HIP point. The CPU shared and independent products were bitwise equal.

Detector-coordinate-corrected flat-pupil comparisons give normalized image
residuals of about 9.5% LLOWFS and 6.9% SCC from 512 to 1024, then 7.6% and
4.2% from 1024 to 2048. These are monotone but not asymptotically converged.
Consequently:

- 512 is the interactive HIL and viewer default, not the precision reference;
- 1024 is the practical higher-fidelity GPU setting;
- 2048 is reserved for selected convergence and calibration checks.

All records use SPIDERSSim `813b4ca`, AdaptiveOpticsProperHIL `8f95ca9`,
AdaptiveOpticsSim `087eb1a`, and Proper `f1aa417`. CUDA ran in detached clean
worktrees. In the local CPU/HIP records, the raw provenance flag for
AdaptiveOpticsSim is dirty only because that checkout contained untracked
review Markdown files; the recorded runtime source revision had no tracked
modifications. SPIDERSSim, AdaptiveOpticsProperHIL, and Proper were clean.
