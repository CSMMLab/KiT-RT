#!/bin/bash
cp -r /mnt /home/rtsn
cd /home/rtsn/
git submodule update --init --recursive
cd code/build/release
rm -rf *
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=True ../../
ninja
cd ../../bin
./unit_tests
