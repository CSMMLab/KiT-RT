# RTSN - an HPC Radio Therapy $`S_n`$ framework
TBD 
## What RTSN is capable of
### Theory
TBD
### Examples
TBD

---

## Build
### Required dependencies
 - Compiler with C++20 support (e.g. gcc/g++ version 9)
 - cmake >= v3.12.4
 - LAPACK
 - OpenMP
 - MPI
 - VTK
 - git
 - ninja or make

### Obtain submodules
Note that an **active internet connection is required for the first build** in order to download the latest versions of the required submodules!
For the first build only, download all submodules:

```bash
git submodule update --init --recursive
```

### Compile the code
**Make** build system (available on most systems)
 
```bash 
cd code/build/release
cmake -DCMAKE_BUILD_TYPE=Release ../../
make 
```
If building in parallel is desired, change the last line to `make -jN`, where `N` optimally is equal to the number of available threads+1.

**Ninja** build system:

```bash 
cd code/build/release
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ../../
ninja
```


The resulting executable will automatically be placed in the `code/bin` folder.

---

## Run
### Local
Execute the compiled binary and hand over a valid *TOML*-styled config file.
Example from inside the `code` directory:

```bash
./bin/RTSN input/example.toml
```

In order to run the code in parallel execute:

```bash
OMP_NUM_THREADS=N mpirun -np J ./bin/RTSN input/example.toml
```

with `N` equal to the number of shared memory threads and `J` equal to the number of distrubuted memory threads.

### BwUniCluster
As VTK is not available on the bwUniCluster, it needs to be installed first. This just needs to be done once. Example:
```bash
module load devel/cmake/3.16
module load compiler/gnu/9.2
wget --no-check-certificate --quiet https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
tar xzf VTK-8.2.0.tar.gz 
mkdir VTK-build
cd VTK-build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCUMENTATION=OFF -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=~/VTK-install ../VTK-8.2.0
make -j10
make install
cd -
rm -r VTK-8.2.0 VTK-build
```

Example for build and run on bwUniCluster:
Get the code
```bash
git clone https://git.scc.kit.edu/rtsn/rtsn.git
cd rtsn/
git submodule init
git submodule update
```
Append `HINTS VTK_INSTALL_DIR` to the `find_package( VTK ... )` line in the CMakeLists.txt. E.g.:
```bash
find_package( VTK REQUIRED COMPONENTS vtkIOGeometry vtkFiltersCore HINTS ~/VTK-install )
```
As a temporary fix, also set the `BLAZE_BLAS_MODE` preprocessor definition to `0`
```bash
add_compile_definitions( BLAZE_BLAS_MODE=0 )
```

Compile it
```bash
module load devel/cmake/3.16
module load compiler/gnu/9.2
module load mpi/openmpi/4.0
cd code/build/release/
cmake -DCMAKE_BUILD_TYPE=Release ../../
make -j10
```

---

## Tests
In order to also compile the unit tests, add `-DBUILD_TESTING=True` to the cmake command, e.g.:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=True ../../
```

The resulting `unit_tests` executable will also be placed in `/bin`.

## Continuous Integration (CI)
Every commit on the master branch will trigger a build test and unit tests.
If either of the tests fail you will get an email, telling you the 'pipeline' has failed. If this happens, visit the 'CI/CD' tab and check what went wrong and try to fix the error as soon as possible.
Creating a merge request from a branch into the master branch will (not tested yet) also trigger the tests and your merge will be rejected in case any test fails.
As a temporary solution, these tests will currently be done on a low power virtual machine running on a 24/7 online PC at Jannick's place (i.e. don't create computationally expensive tests right now as they will fail if exceeding 1GB of memory).
If you add additional libraries to the code, these also need to be added to the test environment, i.e. the respective docker container.
Little guide:
Test the library installation in the docker container
```bash
docker run -it --rm rtsn/test:latest bash
```
Note the steps required and add them to the `Dockerfile` in the `scripts` folder.
Build the new container (takes some time)
```bash
cd docker build -t rtsn/test:latest .
```
or commit your changes to the image (google that procedure).
Push the new image to `hub.docker.com`
```bash
docker push rtsn/test:latest
```
This last step requires a preceeding `docker login`. Ask Jannick for login credentials.

---

## Code structure
![Can't load image](doc/uml.png "UML diagram")
<br/><br/>
The diagram can be modified by opening `doc/uml.zargo` with `ArgoUML`.

---

## Coding style
Please stick to the following coding style for easier code readability:

 - class variables start with an underscore and lowercase letters e.g. `_foo`
 - functions start with a capital letter e.g. `GetSettings()`
 - any variable/function names have capital letters at each individual word e.g. `GetAllCellsAdjacentTo(Cell i)`

Please also use the provided `code/.clang-format` style format to format your code before pushing your latest commits.
Some editors offer to automatically apply the style format upon saving a file (e.g. `Qtcreator`).
