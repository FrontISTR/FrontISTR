FROM ubuntu:18.04 AS build
LABEL maintainer="FrontISTR Commons"

RUN apt-get update \
 && apt-get -y install \
    build-essential    \
    ccache             \
    cmake              \
    curl               \
    gfortran           \
    git                \
    libmetis-dev       \
    libmetis5          \
    libmumps-dev       \
    libopenblas-dev    \
    libopenmpi-dev     \
    libptscotch-dev    \
    ruby               \
    trilinos-all-dev   \
 && apt-get clean      \
 && rm -rf /var/lib/apt/lists/*

FROM build AS document
RUN apt-get update \
 && apt-get install -y \
    doxygen            \
    graphviz           \
    python3            \
    python3-pip        \
 && apt-get clean      \
 && rm -rf /var/lib/apt/lists/*


FROM centos:7 AS cross
RUN yum -y install cmake gcc git gcc-gfortran gcc-c++ binutils.x86_64 glibc-static libstdc++-static epel-release make \\
 && yum -y install mingw-binutils-generic mingw64-cpp mingw64-gcc-c++ mingw64-winpthreads mingw64-gcc mingw64-winpthreads-static mingw64-binutils mingw64-libgomp mingw64-gcc-gfortran mingw64-headers mingw32-cpp mingw32-gcc-c++ mingw32-winpthreads mingw32-gcc mingw32-winpthreads-static mingw32-binutils mingw32-libgomp mingw32-gcc-gfortran mingw32-headers \
 && yum -y clean all       \
 && rm -rf /var/cache/yum/*


