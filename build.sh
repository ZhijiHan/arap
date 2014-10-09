#!/bin/bash
# Change files below for different models.
OFF_FILE_NAME=/home/taodu/research/arap/model/square_21_spikes.off
SELECTION_FILE_NAME=/home/taodu/research/arap/model/square_21_spikes-selection.dmat
DMAT_FILE_NAME=/home/taodu/research/arap/model/square_21_spikes.dmat

cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make clean
make
./demo_bin $OFF_FILE_NAME $SELECTION_FILE_NAME $DMAT_FILE_NAME admm-fixed 500 0.3
