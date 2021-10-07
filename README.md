FrontISTR
======

This is a fork repository of the [FrontISTR](https://github.com/FrontISTR/FrontISTR) tuned for A64FX.
See the [Official home page](https://www.frontistr.com/) for more information about FrontISTR.

# Prerequisites

- FUJITSU Software Compiler Package or Fujitsu Development Studio
- cmake 2.8.0 or later
  - https://cmake.org/download/
- metis 5.1.0
  - http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
- REVOCAP_Refiner 1.1.04
  - https://www.frontistr.com/download/link.php?REVOCAP_Refiner-1.1.04.tar.gz

# Installation
Setup commands and options.
```
LANGDIR="Path to the root directory of Fujitsu compiler"
export PATH=${LANGDIR}/bin:$PATH
export LD_LIBRARY_PATH=${LANGDIR}/lib64:$LD_LIBRARY_PATH
export MPICC=mpifcc
export MPICXX=mpiFCC
export MPIFC=mpifrt
```
Add "px" to MPI[CC|CXX|FC] commands for cross compilation.
```
export MPICC=mpifccpx
export MPICXX=mpiFCCpx
export MPIFC=mpifrtpx
```

## Building metis
```
tar xzvf metis-5.1.0.tar.gz
cd metis-5.1.0
make config prefix="Path to METIS_INSTALL_DIR" cc=${MPICC} openmp=1
make -j 12 install
cd ..
```

## Building REVOCAP_Refiner
Extract archive.
```
tar xzvf REVOCAP_Refiner-1.1.04.tar.gz
cd REVOCAP_Refiner-1.1.04
```

Edit MakefileConfig.in as follows.
```
ARCH = a64fx
CC = ${MPICC}
CFLAGS = -Kcmodel=large -Nlst=t -Kocl -Kfast -Kstrict_aliasing -Krestp=all -Kzfill
CXX = ${MPICXX}
CXXFLAGS = -Kcmodel=large -Nlst=t -Kocl -Kfast -Kstrict_aliasing -Krestp=all -Kzfill
F90 = ${MPIFC}
FFLAGS = -Kcmodel=large -Nlst=t -Kocl -Kfast -Kstrict_aliasing -Krestp=all -Kzfill
AR = ar
ARFLAGS = rsv
LD = ${MPIFC}
LDFLAGS = $(FFLAGS) --linkfortran
LIBPATH =
LIBS = -lm
EXTOBJS =
RM = rm -f
DOXYGEN = doxygen
TAR = tar
```
Make Refiner.
```
make -j 12 Refiner
```

## Building FrontISTR
Clone this repository.
```
$ git clone https://github.com/fujitsu/FrontISTR.git
```

Set path to dependencies and target directory.
```
export METISDIR="Path to metis directory"
export REFINERDIR="Path to REVOCAP_Refiner directory"
export FISTRDIR="Path to FrontISTR install directory"
```

Run cmake and make install.
```
cd FrontISTR
mkdir build
cd build
cmake \
 -DCMAKE_BUILD_TYPE=RELEASE \
 -DCMAKE_C_COMPILER=${MPICC} \
 -DCMAKE_C_FLAGS="-Kcmodel=large -Nlst=t -Kocl -Kfast -Kzfill -Koptmsg=2 -V" \
 -DCMAKE_CXX_COMPILER=${MPICXX} \
 -DCMAKE_CXX_FLAGS="-Kcmodel=large -Nlst=t -Kocl -Kfast -Kzfill -Koptmsg=2 -V" \
 -DCMAKE_Fortran_COMPILER=${MPIFC} \
 -DCMAKE_Fortran_FLAGS="-Kcmodel=large -Nlst=t -Kocl -Kfast -Kzfill -Koptmsg=2 -V" \
 -DWITH_TOOLS=1 \
 -DWITH_MPI=1 \
 -DWITH_OPENMP=1 \
 -DWITH_REFINER=1 \
 -DWITH_METIS=1 \
 -DMETIS_INCLUDE_PATH=${METISDIR}/include \
 -DMETIS_LIBRARIES=${METISDIR}/lib/libmetis.a \
 -DCMAKE_INSTALL_PREFIX=${FISTRDIR} \
 -DREFINER_INCLUDE_PATH=${REFINERDIR}/Refiner \
 -DREFINER_LIBRARIES=${REFINERDIR}/lib/a64fx/libRcapRefiner.a \
 -DCMAKE_Fortran_MODDIR_FLAG=-M \
 -DBLAS_LIBRARIES="-SSL2BLAMP" \
 -DLAPACK_LIBRARIES="-SSL2BLAMP" \
 -DSCALAPACK_LIBRARIES="-SCALAPACK" \
 -DOpenMP_C_FLAGS="-Kopenmp -Kparallel" \
 -DOpenMP_CXX_FLAGS="-Kopenmp -Kparallel" \
 -DOpenMP_Fortran_FLAGS="-Kopenmp -Kparallel" \
 -DCMAKE_EXE_LINKER_FLAGS="--linkfortran -Kopenmp -Kparallel" \
 ..
make -j12 install
```

# Usage
Following are the instructions on how to run the tutorial 02.
```
cd FrontISTR/tutorial/02_elastic_hinge_parallel
${FISTRDIR}/bin/hecmw_part1
mpiexec -n 4 ${FISTRDIR}/bin/fistr1 -t 1
```

# License
See License.txt.
