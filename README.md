# FrontISTR

[![CI Status](https://gitlab.com/frontistr-commons/frontistr/badges/master/pipeline.svg)](https://gitlab.com/frontistr-commons/frontistr/-/pipelines)
[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://manual.frontistr.com/en/)
[![License](https://img.shields.io/badge/license-MIT-green)](License.txt)

**FrontISTR** is an Open-Source Large-Scale Parallel FEM Program for Nonlinear Structural Analysis developed by FrontISTR Commons.

## 📖 Documentation

Detailed manuals are available in multiple languages:
- [English Manual](https://manual.frontistr.com/en/)
- [Japanese Manual](https://manual.frontistr.com/ja/)

## 🚀 Quick Start

### Prerequisites

- C/C++/Fortran compiler (GCC, Intel, etc.)
- CMake (version 3.13.4 or later)
- MPI (optional, for parallel computing)
- OpenMP (optional, for multi-threading)

### Building with CMake (Recommended)

1. **Clone the repository**
   ```bash
   git clone https://gitlab.com/frontistr-commons/frontistr.git
   cd frontistr
   ```

2. **Create build directory**
   ```bash
   mkdir build
   cd build
   ```

3. **Configure and build**
   ```bash
   # Basic build
   cmake ..
   make -j$(nproc)
   
   # With MPI support
   cmake -DWITH_MPI=ON ..
   make -j$(nproc)
   
   # With OpenMP support
   cmake -DWITH_OPENMP=ON ..
   make -j$(nproc)
   
   # Hybrid MPI+OpenMP build
   cmake -DWITH_MPI=ON -DWITH_OPENMP=ON -DWITH_ML=ON -DWITH_MUMPS=ON ..
   make -j$(nproc)
   ```

4. **Run tests**
   ```bash
   ctest --output-on-failure
   ```

### Alternative Build Methods

For users who prefer the traditional approach, makefiles are also available:
```bash
./setup.sh
make
```

## 📁 Directory Structure

```
FrontISTR/
├── CMakeLists.txt          # Main CMake configuration
├── README.md              # This file
├── README.ja.md           # Japanese README
├── LICENSE.txt            # License information
├── VERSION                # Version information
├── build/                 # CMake build directory (created after build)
├── cmake/                 # CMake modules and utilities
├── doc/                   # Documentation sources
├── docker/                # Docker configurations for CI/CD
├── etc/                   # Platform-specific configuration files
├── fistr1/                # FrontISTR main source code
├── hecmw1/                # HEC-MW (mesh and communication library)
├── tests/                 # Test cases and validation
├── tutorial/              # Tutorial examples
└── run_test/              # Test execution directories
```

### Key Components

- **`fistr1/`** - Main FrontISTR finite element solver
- **`hecmw1/`** - HEC-MW library for mesh handling and MPI communication
- **`tests/`** - Comprehensive test suite with various analysis types
- **`tutorial/`** - Step-by-step tutorial examples
- **`cmake/`** - CMake modules for finding dependencies
- **`docker/`** - Container configurations for development and CI

## 🧪 Testing

The project includes extensive test suites covering:
- Serial analysis (`test_serial_*`)
- MPI parallel analysis (`test_mpi_*`)
- OpenMP analysis (`test_openmp_*`)
- Hybrid MPI+OpenMP analysis (`test_hybrid_*`)

Run specific test categories:
```bash
# Run serial tests only
ctest -L serial

# Run MPI tests only
ctest -L mpi

# Run all tests with verbose output
ctest --output-on-failure
```

## 📞 Support

For questions, issues, or contributions:
- **Email**: support@frontistr.com
- **GitLab Issues**: [Report issues](https://gitlab.com/frontistr-commons/frontistr/-/issues)

## 📄 License

Please read [`License.txt`](License.txt) carefully before using this software.

## Acknowledgement

The porting to the GPU was carried out with the support of the Information Technology Center, the University of Tokyo.

---

**Note**: This project uses CMake as the primary build system. While legacy Makefiles are maintained for compatibility, CMake is recommended for new installations and development.
