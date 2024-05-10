FROM ubuntu:24.04 AS build

ARG DEPS="\
        ca-certificates \
        cmake \
        make \
        g++ \
        gemmi \
        gemmi-dev \
        git \
        libboost-program-options-dev \
        libeigen3-dev \
        libfmt-dev \
        libnanoflann-dev \
        libomp-dev \
        libstb-dev \
        nlohmann-json3-dev \
        zlib1g-dev \
        tao-pegtl-dev"

RUN apt-get update && apt-get install -y --no-install-recommends ${DEPS}

ARG REBUILD=foo
RUN git clone --depth 1 https://github.com/sb-ncbr/ChargeFW2.git && \
        cd ChargeFW2 && \
        git checkout master && \
        mkdir build && \
        cd build && \
        cmake .. -DCMAKE_INSTALL_PREFIX=. -DPYTHON_MODULE=OFF && \
        make -j4 && \
        make install

# bundle dependencies
RUN mkdir /dependencies /build
RUN mv /ChargeFW2/build/bin \
        /ChargeFW2/build/lib \
        /ChargeFW2/build/share \
        /build
RUN mv /usr/lib/x86_64-linux-gnu/libgomp.so.1*\
        /usr/lib/x86_64-linux-gnu/libfmt.so* \
        /usr/lib/x86_64-linux-gnu/libboost_program_options.so* \
        /dependencies


FROM ubuntu:24.04 AS app

ENV PATH=/ChargeFW2/build/bin:${PATH}

# copy over the build artifacts
COPY --from=build /build/ /ChargeFW2/build
COPY --from=build /dependencies/* /usr/lib/x86_64-linux-gnu/

# non-root user
ARG UNAME=nonroot
ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} -o ${UNAME} && \
        useradd -m -u ${UID} -g ${GID} -o -s /bin/bash ${UNAME}
ENV HOME=/home/${UNAME}
USER ${UNAME}
WORKDIR ${HOME}

ENTRYPOINT [ "chargefw2" ]
CMD ["--help"]
