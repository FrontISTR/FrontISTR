FROM ubuntu:20.04 AS runtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get -y install libmumps-5.2.1 libgomp1 libtrilinos-ml12 libopenmpi3 libmetis5 openssh-client time \
 && apt-get -y install libmkl-core libmkl-intel-lp64 libmkl-gnu-thread libmkl-scalapack-lp64 libmkl-blacs-openmpi-lp64 \
 && apt-get clean      \
 && rm -rf /var/lib/apt/lists/*

FROM runtime AS build
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
    libmkl-dev         \
 && apt-get clean      \
 && git clone --depth=1  https://gitlab.com/FrontISTR-Commons/REVOCAP_Mesh.git \
 && cd REVOCAP_Mesh \
 && cp config/MakefileConfig.LinuxCluster MakefileConfig.in \
 && make -j4 Refiner \
 && find lib -type f -name "libRcapRefiner*" -exec cp {} /usr/lib/ \; \
 && find . -type f -name "rcapRefiner.h" -exec cp {} /usr/include/ \; \
 && cd ../ && rm -fr ./REVOCAP_Mesh \
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

