FROM ubuntu:25.04 AS build

ARG DEPS="\
        cmake \
        make \
        g++ \
        libboost-program-options-dev \
        libeigen3-dev \
        libfmt-dev \
        libnanoflann-dev \
        libomp-dev \
        nlohmann-json3-dev \
        zlib1g-dev"

RUN apt-get update && apt-get install -y --no-install-recommends ${DEPS}

# Use newer version of Gemmi since Ubuntu currently ships only 0.6.5
ARG GEMMI_VERSION=0.7.3
ADD https://github.com/project-gemmi/gemmi/archive/refs/tags/v${GEMMI_VERSION}.tar.gz .
RUN tar xvzf v0.7.3.tar.gz && \
    cd gemmi-0.7.3 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install

COPY . ChargeFW2
RUN     cd ChargeFW2 && \
        mkdir build && \
        cd build && \
        cmake .. -DCMAKE_INSTALL_PREFIX=. -DPYTHON_MODULE=OFF && \
        make -j$(nproc) && \
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
        /usr/local/lib/libgemmi_cpp.so* \
        /dependencies


FROM ubuntu:25.04 AS app

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
