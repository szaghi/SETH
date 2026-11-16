# Benchmark baselines

Raw outputs of `seth_bench_weno` (16777216 faces, 20 timed passes, the median per variant), three runs per build,
recorded 2026-11-15 at commit `af5f685`, before any kernel added after v0.1.0. `scripts/bench_compare.py` reads a
directory and judges new runs against it: a variant regresses when the median of the new runs exceeds the baseline
median by more than 3 % plus the baseline's run-to-run spread.

| Directory | Build | Machine | Spread of the runs |
|---|---|---|---|
| `rtx4070-nvf-oac/` | `fobis build --mode bench-nvf-oac`, nvfortran 26.9 | NVIDIA RTX 4070 (device 1, `CUDA_VISIBLE_DEVICES=1`), driver 610.88 | ≤ 0.007 ns per face (≤ 0.6 %) |
| `i9-13900kf-core-gnu/` | `fobis build --mode bench-gnu`, gfortran 16.2.0 | one core of an Intel i9-13900KF (`taskset -c 2`), WSL2 | 1 to 17 ns per face (5 to 20 %): a CPU regression of less than that is not detectable here |

```bash
fobis build --mode bench-nvf-oac
for r in 1 2 3; do CUDA_VISIBLE_DEVICES=1 exe-bench/seth_bench_weno > /tmp/run-$r.txt; done
python3 scripts/bench_compare.py notes/bench-baseline/rtx4070-nvf-oac /tmp/run-*.txt
```

A baseline is never edited: a new machine or a deliberate change of a kernel's cost gets a new directory.
