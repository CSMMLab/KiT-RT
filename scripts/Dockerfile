FROM ubuntu:18.04

ENV CXX=/usr/bin/g++ \
    CC=/usr/bin/gcc \
    LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"

RUN apt-get update \
    && apt-get install -qq software-properties-common \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -qq --no-install-recommends \
       gcc-9 \
       g++-9 \
       libopenmpi-dev \
       libblas-dev \
       liblapack-dev \
       git \
       make \
       ninja-build \
       wget \
       libssl-dev \
       libxt-dev \
       libgl1-mesa-dev \
       libgsl-dev \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 100 \
    && update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++-9 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 100 \
    && update-alternatives --install /usr/bin/cc cc /usr/bin/gcc-9 100 \
    && apt-get clean && \
    && apt-get autoremove --purge && \
    && rm -rf /var/lib/apt/lists/* 

RUN wget --no-check-certificate --quiet https://cmake.org/files/v3.16/cmake-3.16.4.tar.gz \
    && tar xzf cmake-3.16.4.tar.gz \
    && cd cmake-3.16.4 \
    && ./bootstrap > /dev/null \
    && make -s -j${nproc} \
    && make -s install > /dev/null \
    && cd - \
    && rm -rf cmake-*
    
RUN wget --no-check-certificate --quiet https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz \
    && tar xzf VTK-8.2.0.tar.gz \
    && mkdir VTK-8.2.0/build \
    && cd VTK-8.2.0/build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCUMENTATION=OFF -DBUILD_TESTING=OFF  ../ \
    && make -s -j${nproc} \
    && make -s install > /dev/null \
    && cd - \
    && rm -rf VTK-*

WORKDIR /home
