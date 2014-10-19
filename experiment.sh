#!/bin/bash
DATA_FOLDER=./data
MODEL_FOLDER=./model
MODEL_NAME=/$1
ITER_NUM=$2
./build/demo_bin $MODEL_FOLDER$MODEL_NAME arap $2 > $DATA_FOLDER$MODEL_NAME-arap-$2.txt

for algorithm in admm-free adapt-admm-free admm-fixed adapt-admm-fixed
  do
    echo $algorithm
    for rho in 0.01 0.02 0.05 0.1 0.2 0.5
      do
        echo $rho
        ./build/demo_bin $MODEL_FOLDER$MODEL_NAME $algorithm $2 $rho > $DATA_FOLDER$MODEL_NAME-$algorithm-$rho.txt
      done
  done
