#!/bin/bash
cp -r /mnt /home/kitrt
cd /home/kitrt/
git submodule update --init --recursive
mkdir build
cd build
rm -rf *
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=True -DBUILD_ML=True ../
ninja
./unit_tests
