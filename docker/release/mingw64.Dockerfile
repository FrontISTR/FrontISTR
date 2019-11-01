
FROM centos:7 AS cross
RUN yum -y install cmake gcc git gcc-gfortran gcc-c++ binutils.x86_64 glibc-static libstdc++-static epel-release make \\
 && yum -y install mingw-binutils-generic mingw64-cpp mingw64-gcc-c++ mingw64-winpthreads mingw64-gcc mingw64-winpthreads-static mingw64-binutils mingw64-libgomp mingw64-gcc-gfortran mingw64-headers mingw32-cpp mingw32-gcc-c++ mingw32-winpthreads mingw32-gcc mingw32-winpthreads-static mingw32-binutils mingw32-libgomp mingw32-gcc-gfortran mingw32-headers \
 && yum -y clean all       \
 && rm -rf /var/cache/yum/*

FROM cross AS mingw64lib
COPY mingw64.cmake /usr/local/mingw64/
RUN export LIB_ROOT=/usr/local/mingw64 \
 && mkdir $LIB_ROOT/include && mkdir $LIB_ROOT/lib \
 && cd /tmp \
 && git clone -b v0.3.7 https://github.com/xianyi/OpenBLAS.git \
 && pushd OpenBLAS \
 && CC=x86_64-w64-mingw32-gcc FC=x86_64-w64-mingw32-gfortran RANLIB=x86_64-w64-mingw32-ranlib HOSTCC=gcc make USE_OPENMP=1 BINARY=64 DYNAMIC_ARCH=1 NO_SHARED=1 -j \
 && make PREFIX=${LIB_ROOT} install \
 && popd \
 && curl -L -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz \
 && tar xvf metis-5.1.0.tar.gz \
 && pushd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h \
 && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir build-static \
 && pushd build-static \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/mingw64.cmake -DOPENMP=ON -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j \
 && make install \
 && popd \
 && popd \
 && curl -L -O http://mumps.enseeiht.fr/MUMPS_5.1.2.tar.gz \
 && tar xvf MUMPS_5.1.2.tar.gz \
 && pushd MUMPS_5.1.2 \
 && cp Make.inc/Makefile.inc.generic.SEQ Makefile.inc \
 && sed -i \
 -e "s|^LAPACK = -llapack|LAPACK = -L${LIB_ROOT}/lib -lopenblas|" \
 -e "s|^LIBBLAS = -lblas|LIBBLAS = -L${LIB_ROOT}/lib -lopenblas|" \
 -e "s|^CC      = cc|CC      = x86_64-w64-mingw32-gcc|"  \
 -e "s|^FC      = f90|FC      = x86_64-w64-mingw32-gfortran|"  \
 -e "s|^FL      = f90|FL      = x86_64-w64-mingw32-gfortran|" \
 -e "s|^OPTF    = -O|OPTF    = -O -fopenmp -DBLR_MT|" \
 -e "s|^OPTC    = -O -I\.|OPTC    = -O -I. -fopenmp|" \
 -e "s|^OPTL    = -O|OPTL    = -O -fopenmp|" Makefile.inc \
 && make RANLIB=x86_64-w64-mingw32-ranlib d -j \
 && cp include/d*.h ${LIB_ROOT}/include \
 && cp lib/*.a ${LIB_ROOT}/lib \
 && cp libseq/*.h ${LIB_ROOT}/include \
 && cp libseq/*.a ${LIB_ROOT}/lib \
 && popd \
 && git clone -b trilinos-release-12-12-1 https://github.com/trilinos/Trilinos.git \
 && pushd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/mingw64.cmake \
  -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_CXX_FLAGS_NONE_OVERRIDE=-fopenmp \
  -DTPL_BLAS_LIBRARIES=$LIB_ROOT/lib/libopenblas.a \
  -DTPL_LAPACK_LIBRARIES=$LIB_ROOT/lib/libopenblas.a \
  -DTPL_METIS_LIBRARIES=$LIB_ROOT/lib/libmetis.a \
  -DTPL_METIS_INCLUDE_DIRS=$LIB_ROOT/include \
  -DTPL_MUMPS_LIBRARIES=$LIB_ROOT/lib/libdmumps.a \
  -DTPL_MUMPS_INCLUDE_DIRS=$LIB_ROOT/include \
  -DTPL_ENABLE_METIS=ON \
  -DTPL_ENABLE_MUMPS=ON \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DTrilinos_ENABLE_Amesos=OFF \
  -DTPL_ENABLE_DLlib=OFF \
  -DMETIS_LIBRARY_DIRS=$LIB_ROOT/lib \
  -DMUMPS_LIBRARY_DIRS=$LIB_ROOT/lib \
  -DBLAS_LIBRARY_DIRS=$LIB_ROOT/lib \
  -DLAPACK_LIBRARY_DIRS=$LIB_ROOT/lib \
  -DBLAS_LIBRARY_NAMES="openblas" \
  -DLAPACK_LIBRARY_NAMES="openblas" \
  -DTrilinos_ENABLE_Fortran=OFF \
  .. \
 && make -j \
 && make install \
 && popd \
 && git clone https://gitlab.com/FrontISTR-Commons/REVOCAP_Mesh.git \
 && cd REVOCAP_Mesh \
 && sed  -i -e "s/ g++/ x86_64-w64-mingw32-g++/" -e "s/AR = ar/AR = x86_64-w64-mingw32-ar/" MakefileConfig.in  \
 && make Refiner -j \
 && find lib -type f -name "libRcapRefiner*" -exec cp {} ${LIB_ROOT}/lib/ \; \
 && find . -type f -name "rcapRefiner.h" -exec cp {} ${LIB_ROOT}/include/ \; \
 && rm -fr /tmp/*



