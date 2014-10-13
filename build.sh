#!/bin/bash
# Change files below for different models.
MODEL_NAME=/home/taodu/research/arap/model/square_21_spikes

cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make clean
make
./demo_bin $MODEL_NAME admm-fixed 500 0.3
