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

