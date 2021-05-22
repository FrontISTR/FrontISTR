FROM debian:buster AS core
RUN apt-get update \
 && apt-get -y install cmake git gfortran make zip unzip curl gcc g++ gfortran mingw-w64-tools gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 gfortran-mingw-w64-x86-64 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

FROM core AS base
COPY toolchain.cmake /usr/x86_64-w64-mingw32/
ENV target=x86_64-w64-mingw32
ENV LIB_ROOT=/usr/x86_64-w64-mingw32
RUN ln -s /usr/x86_64-w64-mingw32 /usr/local/x86_64-w64-mingw32 \
 && git clone --depth 1 https://gitlab.com/FrontISTR-Commons/REVOCAP_Mesh.git && cd REVOCAP_Mesh \
 && cat ./config/MakefileConfig.LinuxCluster|sed  -e "s/ g++/ ${target}-g++/" -e "s/AR = ar/AR = ${target}-ar/" > MakefileConfig.in  \
 && make Refiner -j \
 && find lib -type f -name "libRcapRefiner*" -exec cp {} ${LIB_ROOT}/lib/ \; \
 && find . -type f -name "rcapRefiner.h" -exec cp {} ${LIB_ROOT}/include/ \; \
 && cd .. && rm -fr REVOCAP_Mesh
