# FrontISTR

[![CI Status](https://gitlab.com/frontistr-commons/frontistr/badges/master/pipeline.svg)](https://gitlab.com/frontistr-commons/frontistr/-/pipelines)
[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://manual.frontistr.com/en/)
[![License](https://img.shields.io/badge/license-MIT-green)](License.txt)

FrontISTR is an open-source **large-scale parallel finite element method structural analysis software**. It provides analysis capabilities necessary for industrial and academic practical work, including nonlinear deformation stress, eigenvalue, frequency response, and heat conduction analyses with material nonlinearity, large deformation, and contact. The license is **MIT License**.

* Official Website: [https://www.frontistr.com](https://www.frontistr.com)
* Official Documentation: [https://manual.frontistr.com/en/](https://manual.frontistr.com/en/) (English) / [https://manual.frontistr.com/ja/](https://manual.frontistr.com/ja/) (Japanese)

---

## Key Features

* **Large-scale Parallel**: MPI (distributed memory) + OpenMP (shared memory)
* **Analysis Types**
  * Static analysis, dynamic analysis (implicit/explicit methods)
  * Large deformation, contact, nonlinear materials (elastoplastic, hyperelastic, creep, viscoelastic, etc.)
  * Eigenvalue analysis
  * Frequency response analysis
  * Heat conduction analysis (steady-state/transient)
* **Solvers**: CG / BiCGSTAB / GMRES / GPBiCG / Direct methods (MUMPS / MKL PARDISO)
* **Preconditioners**: SSOR / ILU(0) / AMG (Trilinos-ML), etc.
* **Element Library**: 1st/2nd order solid (TET/PRISM/HEX), plane (triangular/quadrilateral), beam, shell, truss elements, etc.
* **Visualization**: Compatible with ParaView and other tools

---

## System Requirements & Dependencies

* **Compilers**: C/C++/Fortran90 (GCC / Clang / Intel oneAPI, etc.)
* **Required Libraries for Parallel Computing**
  * **MPI** (Open MPI / MPICH, etc.)
  * **METIS** (domain decomposition)
* **Optional Libraries**
  * **BLAS/LAPACK** (used in some features. OpenBLAS/MKL, etc.)
  * **MUMPS** (parallel direct solver) + **ScaLAPACK** (required for MUMPS)
  * **Trilinos-ML** (AMG preconditioning)
  * **Intel MKL** (parallel direct solver with MKL PARDISO)

---

## Build (CMake)

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## Execution

**Single Domain (Thread Parallel)**

```bash
# In the directory where input files are located
fistr1
```

**Domain Decomposition (MPI Parallel)**

```bash
hecmw_part1           # Mesh partitioning
mpirun -n 4 fistr1    # Execute with 4 parallel processes
```

---

## Visualization & Pre/Post Processing

* **ParaView** (VTK file)
* **REVOCAP_PrePost** (mesh generation, preprocessing, execution integration)
* **FreeCAD FEM_FrontISTR** (FrontISTR Workbench): [https://github.com/FrontISTR/FEM_FrontISTR](https://github.com/FrontISTR/FEM_FrontISTR)

---

## Binary Distributions

Download site:

https://www.frontistr.com/download/

---

## Contributing

* Issues and Merge Requests (MR) are accepted on the **[GitLab](https://gitlab.com/FrontISTR-Commons/FrontISTR)** side.
* When proposing changes, please refer to `CONTRIBUTING.md` and various templates (issues, MR).

---

## License

* Main software: **MIT License** (see `License.txt` in the repository)
* Integration libraries (METIS/MUMPS/Trilinos-ML/MKL, etc.) follow their respective licenses.

---

## Acknowledgments

FrontISTR is continuously developed by a development community (FrontISTR Commons) consisting of universities, research institutions, and companies. We appreciate the support from related research and development projects (innovative simulation software, etc.).
The porting to the GPU was carried out with the support of the Information Technology Center, the University of Tokyo.
