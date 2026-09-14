# Tutorial benchmark settings

Edit `cases.json` to select cases and override their input settings.
The same configuration is used for baseline and current runs. Only working
copies are modified; tutorial source files are left unchanged.

- `name`: directory under `tutorial/`.
- `solver`: value after `!SOLVER,METHOD=`. Omit to keep the tutorial setting.
- `solver_parameters`: replacement data lines after `!SOLVER`. Omit to keep
  the existing data lines.
- `requires`: required CMake option, such as `WITH_MKL`.
- `serial` / `mpi`: overrides for non-MPI / MPI **builds**, independent of
  the process count used to run them. OpenMP-only builds use `serial`.
- Top-level `output_type`: visualization output format for all cases.

Solver and output overrides are recorded in the JSON report and reported in
English in the final summary. Notes are generated from these settings rather
than maintained separately in the runner.

## Measurements

For each case, measure baseline then current. If the first successful baseline
measurement takes at most 15 seconds, repeat that pair twice (three runs each).
Compare medians and retain every sample and its logs. Any failed sample makes
that case fail and excludes the pair from timing totals.

Single-run comparisons keep the configured thresholds (by default warning:
15% and 1s; critical: 10% and 5s). Three-run comparisons use warning: 5% and
0.5s; critical: 10% and 1s. Both percentage and absolute limits must be met.
Changes below the thresholds are labeled `WITHIN THRESHOLD`; improvements
use the same magnitude thresholds. Critical regressions fail the job only
when `BENCHMARK_FAIL_ON_REGRESSION` is enabled.

`UNSTABLE` means the sample range is at least 10% of the median and at least
0.5s. This is a separate annotation, not a reason to suppress regression checks.
Cached baselines cannot be interleaved with new measurements; tighter thresholds
apply only when both revisions have three samples.
