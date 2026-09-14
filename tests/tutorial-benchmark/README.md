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
