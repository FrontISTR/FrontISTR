FROM debian:buster AS build
RUN apt-get update \
 && apt-get -y install cmake git gfortran make zip unzip curl \
 && apt-get -y install gcc g++ gfortran \
 && apt-get -y install mingw-w64-tools gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 gfortran-mingw-w64-x86-64 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

FROM build AS lib1
COPY toolchain.cmake /usr/local/x86_64-w64-mingw32/
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && mkdir $LIB_ROOT/include && mkdir $LIB_ROOT/lib \
 && cd /tmp \
 && git clone --depth 1 -b v0.3.13 https://github.com/xianyi/OpenBLAS.git \
 && cd OpenBLAS \
 && CC=${target}-gcc FC=${target}-gfortran RANLIB=${target}-ranlib HOSTCC=gcc LDFLAGS=-fopenmp make USE_OPENMP=1 BINARY=64 DYNAMIC_ARCH=1 NO_SHARED=1 -j \
 && make PREFIX=${LIB_ROOT} install \
 && cd .. \
 && curl -L -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz \
 && tar xvf metis-5.1.0.tar.gz \
 && cd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h \
 && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir build-static \
 && cd build-static \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/toolchain.cmake -DOPENMP=ON -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j \
 && make install \
 && cd .. \
 && cd .. \
 && git clone --depth 1 https://gitlab.com/FrontISTR-Commons/REVOCAP_Mesh.git \
 && cd REVOCAP_Mesh \
 && cat ./config/MakefileConfig.LinuxCluster|sed  -e "s/ g++/ ${target}-g++/" -e "s/AR = ar/AR = ${target}-ar/" > MakefileConfig.in  \
 && make Refiner -j \
 && find lib -type f -name "libRcapRefiner*" -exec cp {} ${LIB_ROOT}/lib/ \; \
 && find . -type f -name "rcapRefiner.h" -exec cp {} ${LIB_ROOT}/include/ \; \
 && rm -fr /tmp/*


FROM lib1 AS lib2
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && curl -L -O http://mumps.enseeiht.fr/MUMPS_5.4.0.tar.gz \
 && tar xvf MUMPS_5.4.0.tar.gz \
 && cd MUMPS_5.4.0 \
 && cp Make.inc/Makefile.inc.generic.SEQ Makefile.inc \
 && make -C src build_mumps_int_def.o build_mumps_int_def \
 && sed -i \
 -e "s|^LAPACK.*$|LAPACK = -L${LIB_ROOT}/lib -lopenblas|" \
 -e "s|^LIBBLAS.*$|LIBBLAS = -L${LIB_ROOT}/lib -lopenblas|" \
 -e "s|^CC.*$|CC = ${target}-gcc|"  \
 -e "s|^FC.*$|FC = ${target}-gfortran|"  \
 -e "s|^FL.*$|FL = ${target}-gfortran|" \
 -e "s|^OPTF.*$|OPTF = -O -fopenmp -DBLR_MT|" \
 -e "s|^OPTC.*$|OPTC = -O -I. -fopenmp|" \
 -e "s|^OPTL.*$|OPTL = -O -fopenmp|" Makefile.inc \
 && make RANLIB=${target}-ranlib d -j \
 && cp include/d*.h ${LIB_ROOT}/include \
 && cp lib/*.a ${LIB_ROOT}/lib \
 && cp libseq/*.h ${LIB_ROOT}/include \
 && cp libseq/*.a ${LIB_ROOT}/lib \
 && cd .. \
 && rm -fr /tmp/*

FROM lib2 AS lib
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && git clone --depth 1 -b trilinos-release-13-0-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake \
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
 && cd .. \
 && rm -fr /tmp/*

FROM lib2 AS lib-trilinos12
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake \
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
 && cd .. \
 && rm -fr /tmp/*

FROM lib2 AS lib-metis4
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && curl -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz \
 && tar xvf metis-4.0.3.tar.gz \
 && cd metis-4.0.3 \
 && sed -e 's/CC = cc/CC = x86_64-w64-mingw32-gcc/' -i Makefile.in \
 && sed -e 's/COPTIONS = /COPTIONS = -D__VC__/' -i Makefile.in \
 && make -C Lib -j \
 && cp libmetis.a $LIB_ROOT/lib/ \
 && find Lib -name "*.h"|xargs -i cp {} $LIB_ROOT/include/ \
 && cd .. \
 && git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake \
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
 && cd .. \
 && rm -fr /tmp/*

FROM lib1 AS mpilib2
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && curl -L -O https://www.frontistr.com/files/MS-MPI_v10.1.2.tar.gz \
 && tar xvf MS-MPI_v10.1.2.tar.gz && cp -r ./msmpi/bin ./msmpi/include ./msmpi/lib ${LIB_ROOT} \
 && cd ${LIB_ROOT}/lib/ && gendef msmpi.dll && x86_64-w64-mingw32-dlltool -d msmpi.def -l libmsmpi.a -D msmpi.dll && rm msmpi.def && cd - \
 && curl -L -O http://www.netlib.org/scalapack/scalapack-2.1.0.tgz \
 && tar xvf scalapack-2.1.0.tgz \
 && cd scalapack-2.1.0 \
 && sed -e "s/mpif90/x86_64-w64-mingw32-gfortran/g" -e "s/mpicc/x86_64-w64-mingw32-gcc/"  -e "s/ranlib/x86_64-w64-mingw32-ranlib/" -e "s/ar/x86_64-w64-mingw32-ar/" -e "s|-O3|-O3 -I${LIB_ROOT}/include|g" SLmake.inc.example > SLmake.inc \
 && make lib -j && cp libscalapack.a ${LIB_ROOT}/lib \
 && cd .. \
 && curl -L -O http://mumps.enseeiht.fr/MUMPS_5.4.0.tar.gz \
 && tar xvf MUMPS_5.4.0.tar.gz \
 && cd MUMPS_5.4.0 \
 && cp Make.inc/Makefile.inc.generic Makefile.inc \
 && make -C src build_mumps_int_def.o build_mumps_int_def \
 && sed \
 -e "s|^LAPACK = -llapack|LAPACK = -L${LIB_ROOT}/lib -lopenblas|" \
 -e "s|^LIBBLAS = -lblas|LIBBLAS = -L${LIB_ROOT}/lib -lopenblas|" \
 -e "s|^CC.*$|CC = ${target}-gcc|"  \
 -e "s|^FC.*$|FC = ${target}-gfortran|"  \
 -e "s|^FL.*$|FL = ${target}-gfortran|" \
 -e "s|^INCPAR.*$|INCPAR = -I${LIB_ROOT}/include|" \
 -e "s|^SCALAP.*$|SCALAP = -lscalapack|" \
 -e "s|^LIBPAR.*$|LIBPAR = \$(SCALAP) \$(LAPACK) ${LIB_ROOT}/lib/libmsmpi.a|" \
 -e "s|^OPTF.*$|OPTF = -O -fopenmp -DBLR_MT|" \
 -e "s|^OPTC.*$|OPTC = -O -I. -fopenmp|" \
 -e "s|^OPTL.*$|OPTL = -O -fopenmp|" -i Makefile.inc \
 && make RANLIB=${target}-ranlib d -j \
 && cp include/d*.h ${LIB_ROOT}/include \
 && cp lib/*.a ${LIB_ROOT}/lib \
 && cd .. \
 && rm -fr /tmp/*


FROM mpilib2 AS mpilib
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && git clone --depth 1 -b trilinos-release-13-0-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake \
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
  -DTPL_ENABLE_MPI=ON \
  .. \
 && make -j \
 && make install \
 && cd .. \
 && rm -fr /tmp/*

FROM mpilib2 AS mpilib-trilinos12
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake \
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
  -DTPL_ENABLE_MPI=ON \
  .. \
 && make -j \
 && make install \
 && cd .. \
 && rm -fr /tmp/*

FROM mpilib2 AS mpilib-metis4
RUN export target=x86_64-w64-mingw32 \
 && export LIB_ROOT=/usr/local/$target \
 && cd /tmp \
 && cd /tmp \
 && curl -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz \
 && tar xvf metis-4.0.3.tar.gz \
 && cd metis-4.0.3 \
 && sed -e 's/CC = cc/CC = x86_64-w64-mingw32-gcc/' -i Makefile.in \
 && sed -e 's/COPTIONS = /COPTIONS = -D__VC__/' -i Makefile.in \
 && make -C Lib -j \
 && cp libmetis.a $LIB_ROOT/lib/ \
 && find Lib -name "*.h"|xargs -i cp {} $LIB_ROOT/include/ \
 && cd .. \
 && git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git \
 && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake \
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
  -DTPL_ENABLE_MPI=ON \
  .. \
 && make -j \
 && make install \
 && cd .. \
 && rm -fr /tmp/*


