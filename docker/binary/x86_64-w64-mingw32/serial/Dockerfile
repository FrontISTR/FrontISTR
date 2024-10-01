FROM registry.gitlab.com/frontistr-commons/frontistr/x86_64-w64-mingw32/base:serial AS lib2
RUN curl -L https://ftp.mcs.anl.gov/pub/petsc/externalpackages/MUMPS_5.4.0.tar.gz | tar zxv && cd MUMPS_5.4.0 \
 && cp Make.inc/Makefile.inc.generic.SEQ Makefile.inc \
 && make -C src build_mumps_int_def.o build_mumps_int_def \
 && sed \
 -e "s|^CC.*$|CC = ${target}-gcc|"  \
 -e "s|^FC.*$|FC = ${target}-gfortran|"  \
 -e "s|^FL.*$|FL = ${target}-gfortran|" \
 -e "s|^INCPAR.*$|INCPAR = -I${LIB_ROOT}/include|" \
 -e "s|^OPTF.*$|OPTF = -O|" \
 -e "s|^OPTC.*$|OPTC = -O -I.|" \
 -e "s|^OPTL.*$|OPTL = -O|" -i Makefile.inc \
 && make RANLIB=${target}-ranlib prerequisites libseqneeded -j \
 && make -C src RANLIB=${target}-ranlib all -j \
 && cp include/*.h ${LIB_ROOT}/include && cp lib/*.a ${LIB_ROOT}/lib && cp libseq/*.h ${LIB_ROOT}/include && cp libseq/*.a ${LIB_ROOT}/lib \
 && cd .. && rm -fr MUMPS_5.4.0

FROM lib2 AS openblas
RUN git clone --depth 1 -b trilinos-release-13-0-1 https://github.com/trilinos/Trilinos.git && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && sed -e '1s/^/#include <windows.h>\n/' -e '1s/^/#include <unistd.h>\n/' -i packages/ml/src/Utils/ml_epetra_utils.cpp \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DTrilinos_ENABLE_OpenMP=OFF \
  -DBUILD_SHARED_LIBS=OFF -DTPL_ENABLE_DLlib=OFF \
  -DTPL_ENABLE_METIS=ON -DTPL_ENABLE_MUMPS=OFF \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DTrilinos_ENABLE_Amesos=ON \  
  -DBLAS_LIBRARY_NAMES="openblas" -DLAPACK_LIBRARY_NAMES="openblas" \
  -DHAVE_GCC_ABI_DEMANGLE=1 -DHAVE_TEUCHOS_BLASFLOAT=1 -DHAVE_TEUCHOS_LAPACKLARND=1 \
  -DTrilinos_ENABLE_Fortran=OFF -DCMAKE_CXX_FLAGS=-fpermissive .. \
 && make -j && make install \
 && cd ../.. && rm -fr Trilinos

FROM lib2 AS openblas-trilinos12
RUN git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && sed -e '1s/^/#include <windows.h>\n/' -e '1s/^/#include <unistd.h>\n/' -i packages/ml/src/Utils/ml_epetra_utils.cpp \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DTrilinos_ENABLE_OpenMP=OFF \
  -DBUILD_SHARED_LIBS=OFF -DTPL_ENABLE_DLlib=OFF \
  -DTPL_ENABLE_METIS=ON -DTPL_ENABLE_MUMPS=OFF \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DTrilinos_ENABLE_Amesos=ON \  
  -DBLAS_LIBRARY_NAMES="openblas" -DLAPACK_LIBRARY_NAMES="openblas" \
  -DHAVE_GCC_ABI_DEMANGLE=1 -DHAVE_TEUCHOS_BLASFLOAT=1 -DHAVE_TEUCHOS_LAPACKLARND=1 \
  -DTrilinos_ENABLE_Fortran=OFF -DCMAKE_CXX_FLAGS=-fpermissive .. \
 && make -j && make install \
 && cd ../.. && rm -fr Trilinos

FROM lib2 AS openblas-metis4
RUN curl -L https://github.com/mfem/tpls/raw/gh-pages/metis-4.0.3.tar.gz | tar zxv && cd metis-4.0.3 \
 && sed -e 's/CC = cc/CC = x86_64-w64-mingw32-gcc/' -i Makefile.in && sed -e 's/COPTIONS = /COPTIONS = -D__VC__/' -i Makefile.in \
 && make -C Lib -j && cp libmetis.a $LIB_ROOT/lib/ \
 && find Lib -name "*.h"|xargs -i cp {} $LIB_ROOT/include/ \
 && cd .. && rm -fr metis-4.0.3

FROM lib2 AS mkl
RUN git clone --depth 1 -b trilinos-release-13-0-1 https://github.com/trilinos/Trilinos.git && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && sed -e '1s/^/#include <windows.h>\n/' -e '1s/^/#include <unistd.h>\n/' -i packages/ml/src/Utils/ml_epetra_utils.cpp \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DTrilinos_ENABLE_OpenMP=OFF \
  -DBUILD_SHARED_LIBS=OFF -DTPL_ENABLE_DLlib=OFF \
  -DTPL_ENABLE_METIS=ON -DTPL_ENABLE_MUMPS=OFF \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DTrilinos_ENABLE_Amesos=ON \  
  -DBLAS_LIBRARY_NAMES="mkl_rt;mkl_intel_lp64_dll;mkl_intel_thread_dll;mkl_core_dll" -DLAPACK_LIBRARY_NAMES="mkl_rt;mkl_intel_lp64_dll;mkl_intel_thread_dll;mkl_core_dll" \
  -DHAVE_GCC_ABI_DEMANGLE=1 -DHAVE_TEUCHOS_BLASFLOAT=1 -DHAVE_TEUCHOS_LAPACKLARND=1 \
  -DTrilinos_ENABLE_Fortran=OFF -DCMAKE_CXX_FLAGS=-fpermissive .. \
 && make -j && make install \
 && cd ../.. && rm -fr Trilinos

FROM lib2 AS mkl-trilinos12
RUN git clone --depth 1 -b trilinos-release-12-18-1 https://github.com/trilinos/Trilinos.git && cd Trilinos \
 && sed -i -e "s/git.cmd/git/" ./cmake/tribits/core/package_arch/TribitsConstants.cmake \
 && sed -e '1s/^/#include <windows.h>\n/' -e '1s/^/#include <unistd.h>\n/' -i packages/ml/src/Utils/ml_epetra_utils.cpp \
 && mkdir build; cd build \
 && cmake \
  -DCMAKE_TOOLCHAIN_FILE=${LIB_ROOT}/toolchain.cmake -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
  -DTrilinos_ENABLE_OpenMP=OFF \
  -DBUILD_SHARED_LIBS=OFF -DTPL_ENABLE_DLlib=OFF \
  -DTPL_ENABLE_METIS=ON -DTPL_ENABLE_MUMPS=OFF \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_ENABLE_TriKota=OFF \
  -DTrilinos_ENABLE_ML=ON \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DTrilinos_ENABLE_OpenMP=ON \
  -DTrilinos_ENABLE_Amesos=ON \  
  -DBLAS_LIBRARY_NAMES="mkl_rt;mkl_intel_lp64_dll;mkl_intel_thread_dll;mkl_core_dll" -DLAPACK_LIBRARY_NAMES="mkl_rt;mkl_intel_lp64_dll;mkl_intel_thread_dll;mkl_core_dll" \
  -DHAVE_GCC_ABI_DEMANGLE=1 -DHAVE_TEUCHOS_BLASFLOAT=1 -DHAVE_TEUCHOS_LAPACKLARND=1 \
  -DTrilinos_ENABLE_Fortran=OFF -DCMAKE_CXX_FLAGS=-fpermissive .. \
 && make -j && make install \
 && cd ../.. && rm -fr Trilinos

FROM lib2 AS mkl-metis4
RUN curl -L https://github.com/mfem/tpls/raw/gh-pages/metis-4.0.3.tar.gz | tar zxv && cd metis-4.0.3 \
 && sed -e 's/CC = cc/CC = x86_64-w64-mingw32-gcc/' -i Makefile.in && sed -e 's/COPTIONS = /COPTIONS = -D__VC__/' -i Makefile.in \
 && make -C Lib -j && cp libmetis.a $LIB_ROOT/lib/ \
 && find Lib -name "*.h"|xargs -i cp {} $LIB_ROOT/include/ \
 && cd .. && rm -fr metis-4.0.3