FROM ubuntu:18.04 AS runtime
RUN apt-get update \
 && apt-get -y install libmumps-5.1.2 libgomp1 libtrilinos-ml12 libopenmpi2 libmetis5 openssh-client \
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
 && apt-get clean      \
 && git clone --depth=1  https://gitlab.com/FrontISTR-Commons/REVOCAP_Mesh.git \
 && cd REVOCAP_Mesh \
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

