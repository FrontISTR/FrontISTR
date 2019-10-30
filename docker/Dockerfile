FROM ubuntu:18.04 AS build
LABEL maintainer="yu"

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

# docker build *MUST* run on the top of this project. see Makefile
FROM build AS fistr1-build
COPY . /src
WORKDIR /src
RUN cmake -B/build -H. -DCMAKE_INSTALL_PREFIX=/src/bin \
 && cmake --build /build --target install -- -j $(nproc)

FROM build AS fistr1
COPY --from=fistr1-build /src/bin /usr
