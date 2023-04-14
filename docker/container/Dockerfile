# docker build *MUST* run on the top of this project. see Makefile
FROM registry.gitlab.com/frontistr-commons/frontistr/build:ubuntu2204 AS fistr1-build
COPY . /src
WORKDIR /src
RUN [ $(uname -m) = x86_64 ] \
      && cmake -B/build -H. -DCMAKE_INSTALL_PREFIX=/usr/local -DWITH_MKL=ON -DWITH_MPI=ON \
      || cmake -B/build -H. -DCMAKE_INSTALL_PREFIX=/usr/local -DWITH_MKL=OFF -DWITH_MPI=ON \
 && cmake --build /build --target install -- -j $(nproc) \
 && rm -fr /src

FROM registry.gitlab.com/frontistr-commons/frontistr/runtime:ubuntu2204 AS fistr1
COPY --from=fistr1-build /usr/local/bin/ /usr/bin/
CMD ["fistr1"]
