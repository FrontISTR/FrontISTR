# Small-strain MITC4 J2 validation

## Scope

The current shell path supports element 741, small strains and rotations,
constant isotropic elasticity and J2 perfect plasticity. Use
`!SOLUTION, TYPE=STATIC, NONLINEAR` and `INFINITESIMAL` in both `!ELASTIC`
and `!PLASTIC`; the second plastic parameter (hardening modulus) must be zero.
The `!SECTION, TYPE=SHELL` data line must specify `thickness, 2`. MITC4 uses
two thickness Gauss points per layer; other requested counts are rejected for
this J2 path instead of silently using two. Existing elastic shell paths are
unchanged. Orthotropic layer definitions are not supported.
The global nonlinear flag enables equilibrium iterations. It does not change
the selected material kinematics into a finite-deformation model.

This is not validation of TL/UL plasticity, hardening, thermal loading,
GPU execution, or finite rotations. MPI checks below cover two-rank small-strain
shell element output, not general MPI correctness or scaling. No existing solid
return-mapping algorithm is changed by this shell extension.

## 4. Thickness strain and output

The local material iteration determines the thickness strain increment so that
the local normal stress is zero. This constitutive thickness strain is retained
in the full strain tensor and transformed for output. It must not be replaced
by the thickness component of the displacement-derived strain.

`check_j2_material.f90` checks elastic uniaxial plane stress, plastic
incompressibility, elastic unloading, and consistent tangents. `run.py` checks
element and nodal results after loading/unloading, then rotates the entire
membrane model and projects the output tensors onto its normal. The global
ZZ stress of an inclined shell need not be zero; the projected normal stress
must be zero. Small-strain output uses engineering shear strains.

The local plane-stress residual tolerance is relative to the larger of the
yield stress and the current maximum stress component, not to a fixed
dimensioned number. Material checks rescale E, yield stress and stress history
by factors from 1e-12 to 1e12 and compare normalized stress, tangent and work,
as well as strain and plastic history. Zero load, elastic loading, plastic
loading and unloading are covered. Solver-level membrane runs additionally
compare nodal, element, layer and material-point output after unit conversion.

`PL_ISTRAIN, ON` writes each shell material-point history as
`PLASTIC_GaussSTRAINk`, without thickness averaging. The one-based index is
`k = ((surface_point-1)*nlayer + layer-1)*nthick + thickness_point`.
MITC4 therefore has 8 values for one layer and 16 for two layers. The output
count is the maximum over all elements and MPI ranks; unused entries are zero.
Elements without shell-layer histories retain their existing Gauss-point data.

`run.py` checks initial, loaded and unloaded membrane output, bending output,
and a model containing both one-layer and two-layer elements. The two layers
have thickness fractions 1/4 and 3/4. Integration-point values must reproduce
the layer +/- averages and the history-point mean equivalent plastic strain;
unused output entries must stay zero. The two stored J2 regression cases also
enable `PL_ISTRAIN` so their integration-point output is compared by CTest.

## 5. Existing output averaging

The J2 extension reads the converged thickness histories and reuses the existing
`get_shell_layer_gauss_average` convention for element strain and stress:
the existing thickness quadrature weights, without adding surface quadrature
weights, layer fractions or Jacobians to that average. Equivalent plastic
strain is the arithmetic mean of the stored material points. Layer +/- element
values retain the arithmetic surface-point mean. They use the outermost
thickness integration points, NOT extrapolated top/bottom material surfaces.
Nodal recovery uses surface extrapolation and the existing adjacent-element
averaging. The mixed-layer nodal averaging issue remains separate.

These representative values must NOT be described as physical-volume or
physical-area averages on distorted elements or with unequal layer thicknesses.
Changing that output definition for both J2 and existing shells is deferred to
a separate Issue, including the earlier `check_j2_weights.f90` trapezoid test.
Raw `ISTRAIN`, `ISTRESS` and `PL_ISTRAIN` values remain available without averaging.

`check_j2_output.f90` instead checks history access, layer values, nodal recovery
and the retained averaging convention using prescribed nonuniform fields at
two unequal layers. It deliberately tests the point mean, not physical volume.
The element stiffness, internal force and potential still use their original
physical integration weights; this separation only concerns result averaging.

## 6. Transverse-shear model

The current implementation retains the work-conjugate transformation

```text
k = 5/6
A = diag(1, 1, 1, 1, sqrt(k), sqrt(k))
strain_work = A * strain_shell
stress_shell = A * stress_work
```

J2 is applied to `stress_work`. After enforcing `stress_work(33)=0`, the
condensed tangent is transformed as `A * D_condensed * A`. Thus the elastic
transverse-shear tangent is kG, while pure transverse-shear yield in the returned
stress is `sqrt(k) * sigma_y / sqrt(3)`. Pure in-plane shear yields at
`sigma_y / sqrt(3)`.

This is a specification of the current corrected shell material model, NOT
proof that the usual physical-stress J2 criterion remains unchanged. The
standard output Mises scalar is computed from the returned stress tensor and
is not the actual yield indicator when transverse shear is present. Changing
to physical-stress J2 with a different treatment of shear correction requires
a separate constitutive decision and further validation. A 5/6 factor in
linear shell theory alone does not justify a nonlinear yield-surface choice.

The material test checks both transverse-shear components, in-plane shear,
elastic slopes, pure-shear yield values, loading work, and finite-difference
tangents for elastic, mixed plastic and unloading states. These verify the
implemented model, not its accuracy against a 3D plasticity solution.

## 7. Unsupported inputs and failures

The setup rejects missing nonlinear equilibrium iterations, finite-strain
material settings, non-perfect-plastic models, orthotropic layer definitions,
unsupported thickness integration counts, temperature-dependent elasticity,
and temperature loads in models containing these J2 shells. This thermal
restriction is model-wide; mixed-model thermal coupling is not implemented.
The material API also checks isotropy, the elastic table, finite positive E and
yield stress, -1 < nu < 1/2, history sizes, and non-finite input/output values.

Local errors distinguish unsupported material, invalid properties/history,
singular thickness tangent, local nonconvergence and non-finite values. The
element propagates the corresponding message instead of reporting every
failure as an unsupported material. These local errors still abort the solve;
automatic cutback/retry for local material failure is NOT implemented.
Unsupported-input and malformed-state rejection paths are exercised. The
finite-strain checks use the actual `KIRCHHOFF` input flag for TL and omit
the kinematics flags for the default UL model. Orthotropy is checked in every
layer, including a material test with an isotropic first layer and orthotropic
second layer. Thickness counts 1, 3 and 5 are rejected; the two-point cases
exercise one-layer and mixed one-/two-layer configurations.
The singular-tangent and local iteration-limit branches have diagnostic messages
but do not yet have a dedicated reproducible failure regression.

## Reproduce

For the current GNU/Linux serial build with LAPACK/BLAS (including WSL), run
from the repository root:

```sh
cmake --build build-shell-ep -j 4
python3 tests/lib/static_lib_shell_mitc4/validation/run.py --build build-shell-ep
ctest --test-dir build-shell-ep -j 4 --output-on-failure
```

The script uses only the Python standard library. It compiles the four Fortran
checks against the supplied build, executes isolated solver inputs and writes
logs plus `summary.json` to a new `build-shell-ep/j2-validation/run-*` directory.
It never replaces checked-in result files. This helper currently assumes GNU
Fortran and a serial LAPACK/BLAS build.

For the optional MPI output comparison, use an MPI-enabled build with the
partitioner and result merger, and the directory printed by `run.py`:

```sh
python3 tests/lib/static_lib_shell_mitc4/validation/check_mpi_output.py \
  --build build-shell-ep-mpi --serial-run build-shell-ep/j2-validation/run-8kn1a756
```

Replace the serial run directory for a new run. This helper uses `mpiexec`,
two ranks and RCB partitioning (no METIS required), and compares global element
IDs after `rmerge`. It writes isolated logs and does not overwrite baselines.

## Local results (2026-10-05)

- All independent checks in `run.py` passed.
- Material FD relative tangent errors: elastic 5.83e-13, plastic 7.30e-10,
  elastic unloading 1.97e-11.
- Unloaded membrane local thickness strain: -0.007947272279808435;
  rotated model: -0.007947272279808437.
- Unloaded membrane nodal/element stress gap: 4.55e-13.
- Bending test: 12 increments completed, maximum 6 Newton iterations;
  support reaction balances the total applied load 0.17.
- Existing serial CTest suite: 103/103 passed. No regression baselines were
  regenerated for the input-validation changes described here.

The membrane test prescribes all nodal DOFs and checks constitutive/output
behavior; it is not evidence of global Newton convergence. The bending test
includes free DOFs and checks nonlinear equilibrium. Stored `.res` files are
computed regression baselines; the independent checks above supply separate
analytical identities and equilibrium checks.

## Integration-point output fix (2026-10-06)

- Single-layer membrane output contains all 8 plastic-strain histories;
  the final value at each point is 0.00794382321525033, not zero.
- A mixed one-/two-layer bending model outputs 16 fields. Layer and history-point
  means agree with the retained convention; entries 9-16 of the one-layer element
  are zero. The two-layer element has nonuniform plastic strain (0 to
  0.002332509933806516), so the check is not limited to a uniform field.
- Only the five J2 regression result files were regenerated to include the
  new output fields. Existing fields changed by at most 2.28e-13.
- All independent checks and 103/103 serial CTest cases passed. MPI execution
  was not tested. The transverse-shear constitutive model was not changed.

## Input and unit-consistency fixes (2026-10-06)

- The new unit-scaling material test failed against the pre-fix library and
  passed after replacing the fixed stress scale in the plane-stress tolerance.
  Across factors 1e-12, 1e-6, 1e6 and 1e12, the largest normalized material
  stress error was 6.26e-16; the solver-level membrane output error was 1.40e-15.
- Orthotropic shell-layer input is rejected even with an isotropic global
  material definition. A two-layer material check also rejects orthotropy
  present only in the second layer.
- Non-two-point thickness inputs are rejected only for the new elastoplastic
  MITC4 path. Both J2 regression meshes now explicitly request two points.
  Variable thickness quadrature has not been implemented.
- Both valid TL (`KIRCHHOFF`) and default UL inputs are rejected with the
  small-strain-only diagnostic. `TOTALLAG` is not used as an input keyword.
- Independent validation passed in `build-shell-ep/j2-validation/run-jk43fap5`.
  All 103 serial CTest cases passed (41.14 s); no reference results were
  regenerated for these fixes. MPI execution was not tested.
- The solid return mapping and the transverse-shear constitutive model remain
  unchanged.

## Quasi-Newton potential connection (2026-10-06)

For the new small-strain J2 shell path, the internal-force quadrature now sums
the layer potential densities into the existing surface-Gauss `strain_energy`
values read by `fstr_get_potential`. The sum includes the surface/thickness
quadrature weights, reference Jacobian, layer fraction, and drilling energy
`alpha * Cv_disp**2 / 2`. It is reset on every trial evaluation. Layer history
keeps the unweighted density; no new persistent arrays or solver branches are
introduced. Existing elastic shell and solid energy paths are unchanged.

`check_j2_potential.f90` checks the energy gradient against all 24 components
of the internal force on a trapezoid with unequal layers. It covers elastic
loading, plastic loading and unloading from committed history, nonzero drilling,
and repeated trial evaluation. Perturbation sizes 1e-7 and 1e-8 give maximum
component-scaled errors of 1.24e-8. The initial elastic potential also agrees
with half the internal-force/displacement product.

The unchanged bending input, with only `METHOD=QUASINEWTON` added, previously
failed at the first increment after 4000 iterations. With the connection it
completes all 12 increments; maximum iteration counts are 2861 (Quasi-Newton)
and 6 (Newton). The maximum displacement difference is 2.60e-6, about 0.0111%
of the maximum Newton displacement. `run.py` compares displacements, rotations
and element/layer/plastic-strain output using a field-scaled tolerance of
5e-4 plus an absolute floor of 1e-9. No regression reference was replaced.

This fixes the missing element potential, not general Quasi-Newton performance.
The driver still has very high iteration counts in this test. It also does not
call `fstr_Update_REACTION_SPC`, unlike the Newton driver, and reports zero
support reactions. Reaction output is recorded, but excluded from the result
agreement check for this documented reason. These existing driver issues are
outside this shell material change; this is not a claim of complete or efficient
Quasi-Newton support. MPI and finite-deformation plasticity are not tested here.

Completion also does not imply the same residual accuracy as Newton: the common
convergence check accepts either the displacement correction or the residual.
The final Quasi-Newton increment in this check meets the correction criterion
while force and moment residuals remain above the requested 1e-8 tolerance.
Iteration-count and reaction-output improvements are deferred to a separate
Issue; the common solver algorithm is not changed here.

## Integration-point strain and stress (2026-10-06)

`ISTRAIN` and `ISTRESS` now read the stored shell thickness-point `strain_out`
and `stress_out`. Each output has six global tensor components, with engineering
shear components for strain. Numbering matches the plastic-strain output:
surface point, layer, then thickness point. A one-layer MITC4 shell has eight
points; a two-layer shell has sixteen. Unused shell point entries are zero-filled
per element, not by clearing every element's output buffer. This applies to
strain, stress and equivalent plastic strain. Existing solid slots retain their
previous behavior, including repeated values beyond an element's own Gauss
count; fixing that legacy padding behavior is outside this change. New fields
beyond the previous model-wide Gauss count have zero entries for solids.

In a six-DOF model, small-strain elastic MITC shells without histories use
`ElementStress_Shell_MITC` at surface Gauss points and thickness quadrature
locations. This also works in purely elastic models when the output is
requested; it does not depend on having a neighboring J2 element.
No material history is allocated or updated for output. The maximum tensor
output count includes the elastic shells' layers. Plastic-strain output retains
the stored-history count and is not expanded by elastic-only layer counts.

`fstr_shell_output_point_count` supplies the supported point count to setup
validation and output. `fstr_get_shell_gauss_output` reads or evaluates one
point's strain and stress together. Its elastic call selects just one surface
point with `surface_gauss_index`; it does not recompute every surface point
for each output slot. The writer only iterates, buffers and registers fields.
The temporary buffer holds one point for all elements and one or two requested
fields, at most 12 doubles per element, and is released after output. There are
no new persistent histories or cached output arrays.

If a six-DOF model contains an element with neither stored thickness histories
nor the supported elastic evaluation path, explicitly requesting `ISTRAIN` or
`ISTRESS` is rejected during setup, before result files are written. Neither
silently zeroing unsupported tensors nor suppressing all requested fields is
used. Unsupported output does not affect runs which do not request these fields.

`run.py` enables both flags in isolated copies of the membrane and bending
inputs. It checks field counts, six-component values, agreement with layer and
element averages, loading/unloading, rotated geometry, stress-unit scaling, and
zero padding in mixed one-/two-layer models. No regression baseline is replaced
for this output change.
Mixed elastic/J2 bending checks now cover output ON/OFF (all existing fields
unchanged), the elastic through-thickness bending distribution, two elastic
layers, and rigidly rotated geometry. The previously zero elastic
`GaussSTRESS1` xx component is -6.516437593340817. The thickness-point mean agrees
with the surface bending stress scaled by 1/sqrt(3); individual surface Gauss
points need not equal the element average. A fully prescribed tetrahedron/brick
plasticity check confirms that all three solid Gauss output families retain
their previous padding behavior with nonzero plastic strain. An unsupported
611 mixture checks setup rejection and the absence of result files.
Strain-only, stress-only and purely elastic cases check field selection and
unchanged existing outputs. `check_shell_output.f90` compares single-point
evaluation with all-point evaluation for elastic MITC3/4/9 shells, rotated
skew geometry and unequal layers; the maximum difference was zero.

Independent validation passed in `build-shell-ep/j2-validation/run-8kn1a756`;
all 103 serial CTest cases passed (42.90 s). No regression baseline was changed
for this output refactoring.

The local GNU/OpenMPI two-rank checks passed in
`build-shell-ep-mpi/j2-validation/mpi-p4vvwxqs`. Eight cases cover plastic bending,
mixed layer counts, elastic/plastic mixtures, elastic multilayers, separate
strain/stress selection and purely elastic output. All element fields are
compared at initial and final states, including unloading where available.
The largest difference scaled by `max(1, abs(reference field))` was 4.06e-10
(purely elastic model); the mixed plastic cases were below 8.38e-14.
A separate pair of disconnected patches puts one-layer and two-layer shells
on different ranks: both ranks write 16 fields, the unused single-layer entries
are zero, and the merged output matches serial exactly. This specifically
checks the global output count, not only a partition with identical local data.

The existing nodal averaging problem when adjacent shells have different layer
counts is deferred to a separate Issue. Its fix is not included here; the new
mixed-layer output checks concern integration-point and element values, not a
claim that mixed-layer nodal averaging is fixed. These deferred Issues have not
been posted to GitLab by this work.

## Result-averaging split (2026-10-06)

The new layer-fraction, surface-quadrature and Jacobian weighting of result
averages is deferred in its entirety. `get_shell_layer_gauss_average` matches
HEAD again; the J2 output recovery uses the existing averaging convention.
The dedicated `check_j2_weights.f90` and its geometry helper are removed from
this change. `check_j2_output.f90` instead checks stored-history output and
nodal recovery using the existing averages. This does not validate physical
volume averages for unequal layers or distorted elements.

The removed implementation and test are preserved outside the repository in
`../FrontISTR_shell_output_average_handoff_20261006`. This is a candidate for
a separate Issue, not a complete fix for all shell output averaging paths.

Comparison with the pre-split `run-2vus0611` found identical displacements,
rotations, reactions and available integration-point values in all 19 completed
validation cases. Single-layer membrane and bending averaged results changed
only at roundoff level (scaled difference below 1.34e-15). The unequal-layer
bending averages intentionally return to the existing convention. All input
and regression reference files are unchanged. Stiffness, internal-force and
potential integration are unchanged, including their physical integration
weights. J2 history output and `ISTRAIN`/`ISTRESS` remain in this change.

## Additional audit (2026-10-06)

The J2 return mapping reads the common `!ELASTIC` properties. Layer-specific
Young moduli and Poisson ratios are not supported by this implementation.
Setup now rejects a shell layer whose constants differ from that table,
including a single-layer mesh/control mismatch. Previously a two-layer input
with Young moduli 40000 and 80000 completed without reporting this mismatch.
The validation now checks both Young-modulus and Poisson-ratio mismatches.
This restriction applies to J2 shells, not to the elastic multilayer path.

`check_j2_potential.f90` now checks all 24 element tangent columns against
internal-force finite differences, tangent symmetry, repeated trial updates,
and a deep history copy/restore without pointer aliasing. Elastic loading,
plastic loading and unloading are covered. At perturbations 1e-7 and 1e-8,
the maximum relative element tangent error was 2.71e-10 in the optimized build
and 2.83e-10 with bounds checks. These are sampled-state checks, not a proof
of convergence for arbitrary load paths.

`run.py` also limits bending to four global Newton iterations and enables
automatic increments. Two failed attempts are cut back successfully; the
analysis finishes in 17 accepted increments with support reaction 0.17.
This exercises global nonconvergence and history restoration, not recovery
from a local constitutive-update error.

The standalone `check_j2_restart.py` compares uninterrupted load/unload runs
with version-6 restarts after loading, for both membrane and bending models.
All final nodal, element and integration-point fields were identical in both
the optimized and bounds-checked builds. Bending unloads over 40 increments.
Running with `--bend-unload-substeps 2` instead reproduces local material
nonconvergence during the first unloading increment. That path aborts rather
than requesting automatic cutback and remains a follow-up limitation; passing
the finer-increment restart check does not resolve it. Quasi-Newton iteration
counts and reaction output also remain separate follow-up items.

Reproduction commands (from the repository root under WSL):

```sh
python3 tests/lib/static_lib_shell_mitc4/validation/run.py --build build-shell-ep
python3 tests/lib/static_lib_shell_mitc4/validation/run.py --build build-shell-ep-check
python3 tests/lib/static_lib_shell_mitc4/validation/check_j2_restart.py --build build-shell-ep
python3 tests/lib/static_lib_shell_mitc4/validation/check_j2_restart.py --build build-shell-ep-check
python3 tests/lib/static_lib_shell_mitc4/validation/check_mpi_output.py --build build-shell-ep-mpi --serial-run build-shell-ep/j2-validation/run-g9wf9cnk
ctest --test-dir build-shell-ep --output-on-failure
```

The check build uses GNU Fortran 13.3, `-O0 -g -fcheck=all -fbacktrace` with
the RELEASE CMake configuration. The DEBUG configuration was not used because
an existing `DEBUG` preprocessor macro conflicts with a parameter of the same
name in `hecmw_precond_SSOR_11.F90`; that unrelated source was not modified.
The optimized build uses `-O3`. Final audit records are:

- Optimized validation: `build-shell-ep/j2-validation/run-g9wf9cnk`.
- Bounds-checked validation: `build-shell-ep-check/j2-validation/run-dggl_utj`.
- Restart: `build-shell-ep/j2-validation/restart-l6fzbpz6` and
  `build-shell-ep-check/j2-validation/restart-lnc_jo4y`.
- Two-rank MPI, eight cases: `build-shell-ep-mpi/j2-validation/mpi-irlw2gw3`.
  Maximum scaled element-field error was 4.06e-10 for the elastic case;
  plastic/mixed cases were below 8.38e-14.
- Final serial CTest: 103/103 passed (156.79 s); the detailed log is
  `build-shell-ep/Testing/Temporary/LastTest.log`. This is a regression run,
  not a controlled performance comparison. The nine J2 baseline input and
  reference files are byte-identical to the pre-audit snapshot.

### Arrays and allocation

Three narrow cleanups avoid unnecessary work without caching more histories:

- Remove the unused inverse-transform argument from the tangent conversion.
- Use assumed-shape output arrays in `NodalStress_ShellJ2`, avoiding copies
  for the caller's noncontiguous work-array sections. The new routine's
  strain/stress temporary warnings are absent in the check build. Existing
  translation-section temporaries in common shell routines remain unchanged.
- Use a fixed two-entry local material status buffer, matching the current
  perfect-plasticity allocation in `fstr_init_gauss`. The optimized
  `ShellMITC_UpdateStress` previously contained direct `malloc/free` calls;
  its final object code contains neither. Extending the hardening scope will
  require revisiting this assumption. No end-to-end speedup is claimed.

On this compiler, `storage_size` plus status payload gives 572 bytes per
Gauss state, excluding allocator overhead. A one-layer J2 MITC4 element keeps
four surface states and eight thickness states, approximately 6864 bytes per
element. The added thickness states account for 4576 bytes of that total.
Thus 100000 elements need about 686 MB for these states alone; a cutback
backup can approximately double that component. This is not total process
memory and varies with compiler layout. No per-element tangent cache was added.

The output writer's temporary tensor buffer is at most 96 bytes per element,
but this is NOT the total result-output memory: `HECMW_result_io_add` copies
and retains each registered field. Eight points with strain, stress and
plastic strain occupy another 832 bytes per element, before other fields,
metadata and allocator overhead. Additional layers/global MPI point-count
padding increase that amount. Large-mesh peak-memory and timing measurements,
heap-leak tools, OpenMP and GPU checks have not been performed in this audit.
