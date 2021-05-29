# docker build *MUST* run on the top of this project. see Makefile
FROM registry.gitlab.com/frontistr-commons/frontistr/build:latest AS fistr1-build
COPY . /src
WORKDIR /src
RUN cmake -B/build -H. -DCMAKE_INSTALL_PREFIX=/src/bin -DWITH_MKL=ON -DWITH_MPI=ON \
 && cmake --build /build --target install -- -j $(nproc)

FROM registry.gitlab.com/frontistr-commons/frontistr/runtime:latest AS fistr1
COPY --from=fistr1-build /src/bin /usr
