FROM ubuntu:22.04

WORKDIR /dependencies

RUN apt-get update && \
    apt-get install -y git cmake libnanoflann-dev \
        libomp-dev libeigen3-dev libfmt-dev \
        libboost-system-dev libboost-program-options-dev libboost-filesystem-dev \
        libfmt-dev nlohmann-json3-dev python3-pybind11 \
        python3-dev g++ gemmi-dev tao-pegtl-dev

RUN git clone --depth 1 https://github.com/krab1k/ChargeFW2.git && \
    cd ChargeFW2 && git checkout master && mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make && make install

ENV PYTHONPATH=/usr/local/lib

WORKDIR /root

ENTRYPOINT [ "bash" ]
