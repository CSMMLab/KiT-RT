[![CI](https://github.com/CSMMLab/KiT-RT/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/CSMMLab/KiT-RT/actions/workflows/c-cpp.yml)
[![Coverage Status](https://coveralls.io/repos/github/CSMMLab/KiT-RT/badge.svg?branch=master)](https://coveralls.io/github/CSMMLab/KiT-RT?branch=master)
[![Documentation Status](https://readthedocs.org/projects/kit-rt/badge/?version=latest)](https://kit-rt.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# KiT-RT - an HPC Radio Therapy simulation framework
The KiT-RT (Kinetic Transport Solver for Radiation Therapy) framework is a high-performance open source platform for radiation transport. Its main focus is on radiotherapy planning in cancer treatment. To enable problem-specific method selection, the framework provides different deterministic solver types. This not only facilitates treatment planning, but also provides tools to investigate various research questions in the field of radiative transfer. This goal is supported by an easily extendable code structure that allows for straightforward implementation of additional methods and techniques.

The documentation can be found [here](https://kit-rt.readthedocs.io/en/develop/index.html). 
The preprint for the KiT-RT journal paper is available [here](https://arxiv.org/abs/2205.08417).
## Contents

* [Capability](#what-kit-rt-is-capable-of)
* [Build](#build)
* [Run](#run)
* [Unit Tests](#unit-tests)
* [Docker](#docker)
* [Tensorflow](#tensorflow)
* [Citation Info](#Cite)

## What KiT-RT is capable of
### Theory
A short description of kinetic theory can be found [here](https://kit-rt.readthedocs.io/en/develop/physics.html).

## Build
### Required dependencies
 - Compiler with C++17 support
 - cmake >= v3.16
 - LAPACK
 - OpenMP
 - MPI
 - python3
 - VTK
 - git
 
### Tensorflow installation
If you choose to enable the machine learning tools via the BUILD_ML option, you need to install the tensorflow C-API:
```
FILENAME=libtensorflow-cpu-linux-x86_64-2.7.0.tar.gz
wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/libtensorflow/${FILENAME}
tar -C /usr/local -xzf ${FILENAME}
ldconfig /usr/local/lib
```

### Python dependencies
- pydicom
- numpy
- pygmsh version 6.1.1 
```bash
  pip install pygmsh==6.1.1
```
 (note that newer versions are not yet supported)


### Obtain submodules
Note that an **active internet connection is required for the first build** in order to download the suitable versions of the required submodules!
For the first build only, download all submodules:

```bash
git submodule update --init --recursive
```

### Compile the code
In case of the **make** build system (available on most systems) run:
 
```bash 
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j
```

---

## Run
Execute the compiled binary by handing over a [valid config file](https://kit-rt.readthedocs.io/en/latest/configFiles.html), e.g.:

```bash
./KiT-RT ../examples/linesource_SN.cfg
```

In order to run the code in parallel execute:

```bash
OMP_NUM_THREADS=N mpirun -np J ./KiT-RT ../examples/linesource_SN.cfg
```

with `N` equal to the number of shared memory threads and `J` equal to the number of distrubuted memory threads.

---

## Unit Tests
After compiling the framework with:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON ../
make -j
```

Unit test can be run with:
```bash
make test
```

---

## Docker
A preconfigured docker container can also be used to run the code.
By running

```bash
docker run --rm -ti -v $(pwd):/home kitrt/test:latest
```


To use the tensorflow backend for neural network based entropy closures,  please use the tensorflow equippend containter by running
```bash
docker run --rm -ti -v $(pwd):/home kitrt/test_ml:latest
```
from the current folder will be mounted to the docker container and the code can be executed without any of the required dependencies.

Bash scripts are provided in the folder tools/CI to get started with the docker environments. To start an interactive docker environment, execute
```bash
docker run -i -t --rm -v $(pwd)/../..:/mnt kitrt/test:latest /bin/bash
```

## Tensorflow
To use the neural network based entropy closure, one needs to install the tensorflow C++ backend. Just use the following commands:
```
FILENAME=libtensorflow-cpu-linux-x86_64-2.7.0.tar.gz
wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/libtensorflow/${FILENAME}
tar -C /usr/local -xzf ${FILENAME}
ldconfig /usr/local/lib
```
and for a gpu based version (you need supported hardware and gpu drivers, see [here](https://www.tensorflow.org/install/gpu) ):
```
FILENAME=libtensorflow-gpu-linux-x86_64-2.7.0.tar.gz
wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/libtensorflow/${FILENAME}
tar -C /usr/local -xzf ${FILENAME}
ldconfig /usr/local/lib
```
## Cite
```
@article{kitrt2023,
    title={KiT-RT: An extendable framework for radiative transfer and therapy},
    author={Kusch, Jonas and Schotth{\"o}fer, Steffen and Stammer, Pia and Wolters, Jannick and Xiao, Tianbai},
    journal={ACM Transactions on Mathematical Software},
    volume={49},
    number={4},
    pages={1--24},
    year={2023},
    publisher={ACM New York, NY}
}
```
