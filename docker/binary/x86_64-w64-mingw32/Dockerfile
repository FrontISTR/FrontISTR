FROM debian:buster AS env1
RUN apt-get update \
 && apt-get -y install gcc g++ gfortran \
 && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN apt-get update \
 && apt-get -y install make cmake git zip unzip curl \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

FROM env1 AS env
RUN apt-get update \
 && apt-get -y install binutils-mingw-w64-x86-64 mingw-w64-common mingw-w64-tools gcc-mingw-w64-base mingw-w64-x86-64-dev \
 && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN apt-get update \
 && apt-get -y install gcc-mingw-w64-x86-64 \
 && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN apt-get update \
 && apt-get -y install g++-mingw-w64-x86-64 gfortran-mingw-w64-x86-64 \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

FROM env AS core1
COPY toolchain.cmake /usr/x86_64-w64-mingw32/
ENV target=x86_64-w64-mingw32
ENV LIB_ROOT=/usr/x86_64-w64-mingw32
RUN ln -s /usr/x86_64-w64-mingw32 /usr/local/x86_64-w64-mingw32 \
 && git clone --depth 1 https://github.com/FrontISTR/REVOCAP_Refiner.git && cd REVOCAP_Refiner \
 && sed  -e "s/ g++/ ${target}-g++/" -e "s/AR = ar/AR = ${target}-ar/" -i MakefileConfig.in  \
 && make Refiner -j4 \
 && find lib -type f -name "libRcapRefiner*" -exec cp {} ${LIB_ROOT}/lib/ \; \
 && find . -type f -name "rcapRefiner.h" -exec cp {} ${LIB_ROOT}/include/ \; \
 && cd .. && rm -fr REVOCAP_Refiner

FROM core1 AS core
RUN curl -L https://www.frontistr.com/files/ci/mkl_2021.2.0.tar.xz | tar Jxv -C ${LIB_ROOT} --wildcards './lib/mkl_avx*' './lib/mkl_core*'
RUN curl -L https://www.frontistr.com/files/ci/mkl_2021.2.0.tar.xz | tar Jkxv -C ${LIB_ROOT} 2>/dev/null || true

FROM core AS serial
RUN git clone --depth 1 -b v0.3.15 https://github.com/xianyi/OpenBLAS.git && cd OpenBLAS \
 && CC=${target}-gcc FC=${target}-gfortran RANLIB=${target}-ranlib HOSTCC=gcc make USE_OPENMP=0 BINARY=64 DYNAMIC_ARCH=1 NO_SHARED=1 -j \
 && make PREFIX=${LIB_ROOT} install \
 && cd .. && rm -fr OpenBLAS
RUN curl -L https://www.phy.bnl.gov/~bviren/spack-mirror/metis/metis-5.1.0.tar.gz | tar zxv && cd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir buildserial && cd buildserial \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/toolchain.cmake -DOPENMP=OFF -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j && make install \
 && cd ../.. && rm -fr metis-5.1.0

FROM core AS thread
RUN git clone --depth 1 -b v0.3.15 https://github.com/xianyi/OpenBLAS.git && cd OpenBLAS \
 && CC=${target}-gcc FC=${target}-gfortran RANLIB=${target}-ranlib HOSTCC=gcc LDFLAGS=-fopenmp make USE_OPENMP=1 BINARY=64 DYNAMIC_ARCH=1 NO_SHARED=1 -j \
 && make PREFIX=${LIB_ROOT} install \
 && cd .. && rm -fr OpenBLAS
RUN curl -L https://www.phy.bnl.gov/~bviren/spack-mirror/metis/metis-5.1.0.tar.gz | tar zxv && cd metis-5.1.0 \
 && sed -i -e "/#include <sys\/resource.h>/d" ./GKlib/gk_arch.h && sed -i -e "/extern int gk_getopt/d" -e "/longopts/d" ./GKlib/gk_getopt.h \
 && mkdir buildserial && cd buildserial \
 && cmake  -DCMAKE_TOOLCHAIN_FILE=$LIB_ROOT/toolchain.cmake -DOPENMP=ON -DGKRAND=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_VERBOSE_MAKEFILE=1  -DGKLIB_PATH=../GKlib -DCMAKE_INSTALL_PREFIX="${LIB_ROOT} " .. \
 && make -j && make install \
 && cd ../.. && rm -fr metis-5.1.0
