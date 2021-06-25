FROM ubuntu:20.04

WORKDIR /dependencies

ENV TZ=US/Central
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y git build-essential cmake ninja-build vim wget \
        lsb-release software-properties-common libomp-dev \
        libgtest-dev libeigen3-dev libboost-all-dev \
        libfmt-dev nlohmann-json3-dev python3-pybind11 \
        python3-dev python3-pip clang

RUN git clone https://github.com/jlblancoc/nanoflann.git && \
    cd nanoflann && mkdir build && cd build && \
    cmake .. -DCMAKE_CXX_COMPILER=clang++ && \
    make && make install

RUN git clone https://github.com/project-gemmi/gemmi.git && \
    cp -r gemmi/include/gemmi /usr/include

RUN git clone https://github.com/danny305/ChargeFW2.git && \
    cd ChargeFW2 && git checkout master && mkdir build && cd build && \
    cmake .. -GNinja \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_BUILD_TYPE=Release && \
    ninja && ninja install

ENV PYTHONPATH=/usr/local/lib

ENTRYPOINT [ "bash" ]
