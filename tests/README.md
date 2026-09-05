[Japanese Version](./README.ja.md)

Run test
---------

Use [ctest][ctest].
We assume here that FrontISTR project is built by [cmake][cmake]:

```
git clone https://gitlab.com/FrontISTR-Commons/FrontISTR
cd FrontISTR   # we call here as FRONTISTR_HOME
cmake . -Bbuild
cmake --build build/ -j $(nproc)
```

Then `fistr1` and other executables are built in `${FRONTISTR_HOME}/build`.
You can run test on this directory:

```
cd build/
ctest
```

### Labels of tests

Tests are managed by ctest's label. There are 4 labels about parallelization:

| label | OpenMP | MPI |
|:------|:------:|:---:|
|serial | OFF    | OFF |
|openmp | ON     | OFF |
|mpi    | OFF    | ON  |
|hybrid | ON     | ON  |

To execute labeled tests, please run [ctest][ctest] with `-L` (`--label-regex`) flag:

```
ctest -L mpi
```

In addition to these parallelization labels, there are labels for "target".
As described below, [cmake][cmake] seeks tests in this directory,
and put a label `analysis/eigen/exK` to tests for `${FRONTISTR_HOME}/tests/analysis/eigen/exK` for example.
To run tests on this directory, please use this label:

```
ctest -L analysis/eigen/exK
```

Because `-L` flag can select by partial match, 

```
ctest -L analysis
```

will executes all tests in `${FRONTISTR_HOME}/tests/analysis`

### Configure test output

As a general technique of [ctest][ctest],

```
ctest -V
```

displays all output of test processes, and

```
ctest --output-on-failure
```

displays the output of failed tests. See `ctest -h` for detail.

Add test
---------

[cmake][cmake] seeks mesh data (`*.msh`) from `${FRONTISTR_HOME}/tests/`, and then registers it as a test target.
This target compares the result of `fistr1` of current build with reference build,
and tests the difference is enough small.
In order to append a new test, you should

1. Create a directory under `${FRONTISTR_HOME}/tests`
    - Directories created in analysis, lib and solver are always included in the test run.
    - Directories created in with_[mkl|mumps|ml] are included in the test run when the cmake option -DWITH_[MKL|MUMPS|ML] is ON.
    - Directories created in _archive are not included in the test run.
2. Put `*.msh` and `*.cnt` files
3. Generate reference result by `${FRONTISTR_HOME}/tests/create_reference.sh`
4. and Confirm this result is correct by your eye

Be sure that because `create_reference.sh` uses executable binary `${FRONTISTR_HOME}/build/fistr/fistr1` by default, you have to build it in advance.
Be sure that because `create_reference_docker.sh` uses official release image, you need the authority for execution of [Docker][docker].

### Frequency response analysis

Frequency response analysis (`!SOLUTION,TYPE=DYNAMIC` with the second item of the first `!DYNAMIC` data line set to `2`)
superposes the modes obtained by a preceding eigen analysis.
When `<mesh>_eigen.cnt` is placed beside `<mesh>.msh`, `test.sh` and `create_reference.sh` run the eigen analysis first
and hand its outputs over to the frequency response analysis.

| handed over | file name |
|:------------|:----------|
| eigenvalues (log)      | `eigen_0.log` |
| eigenvectors (result)  | `<mesh>_eigen.res` |
| time history output    | `<mesh>_dyna.res` |

So `!EIGENREAD` of `<mesh>.cnt` has to specify `eigen_0.log` as its log file.
Both the frequency sweep result `<mesh>.res.0.*` and the time history result `<mesh>_dyna.res.0.*` are compared with the reference.

The tests under `${FRONTISTR_HOME}/tests/analysis/freq` are separated by the condition each of them varies.
A new test case is added to the directory of the condition it exercises, and the other conditions are kept
as the base configuration so that the effect of the varied condition stays isolated.

| directory  | varied condition | cases |
|:-----------|:-----------------|:------|
| `element`  | element type     | `FQ341` `FQ342` `FQ351` `FQ352` `FQ361` `FQ362` |
| `load`     | `!FLOAD`         | `FQL01` real, `FQL02` imaginary, `FQL03` complex on several groups and DOFs, `FQL04` traction on a surface group |
| `boundary` | `!BOUNDARY`      | `FQB01` cantilever, `FQB02` additional sliding support, `FQB03` both ends fixed |
| `modal`    | damping and mode range | `FQM01` mass proportional, `FQM02` mass and stiffness proportional, `FQM03` `!EIGENREAD` starting from the 2nd mode |
| `output`   | `!OUTPUT_RES`    | `FQO01` nodal strain and principal values, `FQO02` elemental values |

The base configuration is a 1.0 x 0.2 x 0.1 cantilever of 361 elements fixed at `FIX`,
loaded at `LOADP` in the z direction, 5 modes, Rayleigh damping of alpha=0 and beta=1.0E-4,
and a frequency sweep of 5 points up to 250 Hz.

The response is complex, so every item of the result file is written twice, with the
label suffixed by `_real` and `_imag`.

[cmake]: https://cmake.org/cmake/help/latest/manual/cmake.1.html
[ctest]: https://cmake.org/cmake/help/latest/manual/ctest.1.html
[docker]: https://www.docker.com/
