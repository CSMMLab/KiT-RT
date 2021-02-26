.. _implementation:

Implementation
------------------------
*****
Build
*****

Required dependencies
=====================
 - Compiler with C++17 support
 - cmake >= v3.12.4
 - LAPACK
 - OpenMP
 - MPI
 - VTK
 - git
 - ninja or make

Obtain submodules
==================
Note that an **active internet connection is required for the first build** in order to download the latest versions of the required submodules!
For the first build only, download all submodules:

.. code-block:: bash 

        git submodule update --init --recursive

Compile the code
================
**Make** build system (available on most systems)
 
 
.. code-block:: bash 

       cd code/build/release
       cmake -DCMAKE_BUILD_TYPE=Release ../../
       make 

If building in parallel is desired, change the last line to `make -jN`, where `N` optimally is equal to the number of available threads+1.

**Ninja** build system:

.. code-block:: bash
 
      cd code/build/release
      cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ../../
      ninja



The resulting executable will automatically be placed in the `code/bin` folder.

----------------------------------------------------------

**********
Run
**********

Local
===========

Execute the compiled binary from the `bin` folder and hand over a valid *TOML*-styled config file.
Example from inside the `code` directory:

```bash
./KiT-RT ../input/example.cfg
```

In order to run the code in parallel execute:

```bash
OMP_NUM_THREADS=N mpirun -np J ./KiT-RT ../input/example.cfg
```

with `N` equal to the number of shared memory threads and `J` equal to the number of distrubuted memory threads.

BwUniCluster
==============

As VTK is not available on the bwUniCluster, it needs to be installed first. This just needs to be done once. Example:

.. code-block:: bash

           module load devel/cmake/3.16
           module load compiler/gnu/9.2
           wget --no-check-certificate --quiet https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
           tar xzf VTK-8.2.0.tar.gz 
           mkdir VTK-build
           cd VTK-build
           cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCUMENTATION=OFF -DBUILD_TESTING=OFF - 
           DCMAKE_INSTALL_PREFIX=~/VTK-install ../VTK-8.2.0
           make -j
           make install
           cd -
           rm -r VTK-8.2.0 VTK-build


Example for build and run on bwUniCluster:
Get the code

.. code-block:: bash

          git clone https://git.scc.kit.edu/rtsn/rtsn.git KiT-RT
          cd KiT-RT/
          git submodule init
          git submodule update

Append ``HINTS VTK_INSTALL_DIR` to the ``find_package( VTK ... )`` line in the CMakeLists.txt. E.g.:

.. code-block:: bash

          find_package( VTK REQUIRED COMPONENTS vtkIOGeometry vtkFiltersCore HINTS ~/VTK-install )


Compile it

.. code-block:: bash

        module load devel/cmake/3.16
        module load compiler/gnu/9.2
        module load mpi/openmpi/4.0
        cd code/build/release/
        cmake -DCMAKE_BUILD_TYPE=Release ../../
        make -j


---------------------------------------------------------------

Tests
============================

After compiling the framework as described above just run:

.. code-block:: bash
        
		make test


The ``unit_tests`` executable will also be placed in in the build folder.




