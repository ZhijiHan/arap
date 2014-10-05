#!/bin/bash
# Change files below for different models.
OFF_FILE_NAME=/home/taodu/research/arap/model/decimated-knight.off
DMAT_FILE_NAME=/home/taodu/research/arap/model/decimated-knight-selection.dmat

cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make clean
make
./demo_bin $OFF_FILE_NAME $DMAT_FILE_NAME 120 admm 500 0.3
