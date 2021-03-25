# FrontISTR for SX-Aurora TSUBASA

## How to Build

```
make -f ./Makefile.SX
```

- If successful, Directory `fistr_VH` and `hecmw_VE` are created.
- Specify the binary file of `hecmw_VE/hecmw_ve` in the environment variable `FRONTISTR_VEOBJ`.
- The analysis is performed using fistr in `fistr_VH/fistr1`.

## System Requirements
- ncc 3.1.1 or higher
- nfort 3.1.1 or higher
- gfortran 4.8.5 or higher
- [aveo v0.9.15](https://github.com/SX-Aurora/aveo/tree/v0.9.15)

## Restrictions
- MPI is not supported
- Multi-VE is not supported
- Only CG method can be used
- Preprocessing has not been confirmed to work
- Accurate FrontISTR time logs is not supported (There is a lag in synchronizing variables that measure time between VH and VE. If you need accurate time measurement, insert time measurement and time output into code)
- Contents of the output file is not supported when the program terminates abnormally
