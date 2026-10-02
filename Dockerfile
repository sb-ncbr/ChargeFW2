FROM ubuntu:26.04 AS build

ARG DEPS="\
        cmake \
        make \
        g++ \
        libboost-program-options-dev \
        libeigen3-dev \
        libnanoflann-dev \
        libomp-dev \
        gemmi \
        libgemmi-dev \
        nlohmann-json3-dev\
        python3-pybind11"

RUN apt-get update && apt-get install -y --no-install-recommends ${DEPS}

ARG PORTABLE=OFF
COPY . ChargeFW2
RUN     cd ChargeFW2 && \
        mkdir build && \
        cd build && \
        cmake .. \
        -DCMAKE_INSTALL_PREFIX=. \
        -DPYTHON_MODULE=ON \
        -DCMAKE_BUILD_TYPE=Release \
        -DCHARGEFW2_PORTABLE=${PORTABLE} && \
        make -j$(nproc) && \
        make install

# bundle dependencies
RUN mkdir /dependencies /build
RUN mv /ChargeFW2/build/bin \
        /ChargeFW2/build/lib \
        /ChargeFW2/build/share \
        /build
RUN mv /usr/lib/x86_64-linux-gnu/libgomp.so.1*\
        /usr/lib/x86_64-linux-gnu/libboost_program_options.so* \
        /dependencies

FROM ubuntu:26.04 AS app

ENV PATH=/ChargeFW2/bin:${PATH}

# copy over the build artifacts
COPY --from=build /build/ /ChargeFW2/
COPY --from=build /dependencies/* /usr/lib/x86_64-linux-gnu/

# Setup ENV variables for use with Python bindings
ENV CHARGEFW2_INSTALL_DIR=/ChargeFW2/
ENV LD_LIBRARY_PATH=${CHARGEFW2_INSTALL_DIR}/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=${CHARGEFW2_INSTALL_DIR}/lib

USER ubuntu

ENTRYPOINT [ "chargefw2" ]
CMD ["--help"]
