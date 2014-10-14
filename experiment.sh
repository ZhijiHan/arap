#!/bin/bash
DATA_FOLDER=./data
MODEL_FOLDER=./model
MODEL_NAME=/$1
./build/demo_bin $MODEL_FOLDER$MODEL_NAME arap 100 > $DATA_FOLDER$MODEL_NAME-arap-100.txt

for algorithm in admm-free adapt-admm-free admm-fixed adapt-admm-fixed
  do
    echo $algorithm
    for rho in 1 2 5 10 20 50
      do
        echo $rho
        ./build/demo_bin $MODEL_FOLDER$MODEL_NAME $algorithm 100 $rho > $DATA_FOLDER$MODEL_NAME-$algorithm-$rho.txt
      done
  done
