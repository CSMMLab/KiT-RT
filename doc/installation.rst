.. _installation:

Installation
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

---------------------------------------------------------------

Tests
============================

After compiling the framework as described above just run:

.. code-block:: bash
        
		make test


The ``unit_tests`` executable will also be placed in in the build folder.
