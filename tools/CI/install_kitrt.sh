#!/bin/bash
cp -r /mnt /home/kitrt
cd /home/kitrt/
git submodule update --init --recursive
mkdir build
cd build
rm -rf *
cmake -G Ninja -DCMAKE_BUILD_UNITY=True -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=False -DBUILD_ML=True ../
ninja
cd ../examples
../build/KiT-RT checkerboard_MN.cfg
cp -r /home/result /mnt/result
