FROM ubuntu:20.04

ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib" \
    PYTHONPATH=/usr/local/gmsh/lib:$PYTHONPATH \
    PATH=/usr/local/gmsh/bin:$PATH 

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -qq \
       gcc \
       g++ \
       libopenmpi-dev \
       openmpi-bin \
       libblas-dev \
       liblapack-dev \
       git \
       make \
       ninja-build \
       cmake \
       wget \
       ssh \
       libssl-dev \
       libxt-dev \
       libgl1-mesa-dev \
       libglu1 \
       libxrender1 \
       libxcursor-dev \
       libxft-dev \
       libxinerama-dev \
       python3 \
       python3-pip \
       doxygen \
    && apt-get clean \
    && apt-get autoremove --purge \
    && rm -rf /var/lib/apt/lists/* 

RUN cd /usr/local \
    && wget -nc --quiet  http://gmsh.info/bin/Linux/gmsh-4.7.0-Linux64-sdk.tgz \
    && tar xzf gmsh-4.7.0-Linux64-sdk.tgz \
    && mv gmsh-4.7.0-Linux64-sdk gmsh \
    && rm gmsh-4.7.0-Linux64-sdk.tgz 

RUN wget -nc --no-check-certificate --quiet https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz \
    && tar xzf VTK-8.2.0.tar.gz \
    && mkdir VTK-8.2.0/build \
    && cd VTK-8.2.0/build \
    && cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCUMENTATION=OFF -DBUILD_TESTING=OFF  ../ \
    && ninja \
    && ninja install > /dev/null \
    && cd - \
    && rm -rf VTK-*

RUN pip3 install numpy pygmsh==6.1.1 Pillow pydicom gcovr sphinx_rtd_theme breathe

RUN echo "Installation successful"
WORKDIR /home