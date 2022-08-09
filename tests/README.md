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
    - Directories created in _with_** are included in the test run when the cmake option -DWITH_** is ON.
    - Directories created in _archive are not included in the test run.
2. Put `*.msh` and `*.cnt` files
3. Generate reference result by `${FRONTISTR_HOME}/tests/create_reference.sh`
4. and Confirm this result is correct by your eye

Be sure that because `create_reference.sh` uses executable binary `${FRONTISTR_HOME}/build/fistr/fistr1` by default, you have to build it in advance.
Be sure that because `create_reference_docker.sh` uses official release image, you need the authority for execution of [Docker][docker].

[cmake]: https://cmake.org/cmake/help/latest/manual/cmake.1.html
[ctest]: https://cmake.org/cmake/help/latest/manual/ctest.1.html
[docker]: https://www.docker.com/
