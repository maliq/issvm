#!/bin/bash
#dataset_dir=/home/ricky/Cotter-Code
#calculate_error=/home/ricky/Cotter-Code/issvm/computeMissRate.py
dataset_dir=..
calculate_error=${dataset_dir}/issvm/computeMissRate.py

METHOD=sparsifier
DATASET=w8a
ITERATIONS=100000
TRAIN_DATA=Web/w8a.original.train.01.txt
TEST_DATA=Web/w8a.original.test.01.txt
VAL_DATA=Web/w8a.original.test.01.txt.val.
TEST_DATA=Web/w8a.original.test.01.txt.test.
INIT_DIR=$METHOD/init
MODEL_DIR=$METHOD/model
TEST_DIR=$METHOD/test
K=1.0
#declare -a norms=(12.36 18.57 27.20 45.99 100.50 267.46 762.54 1895.5 3106.05 3142.15)
norms_init=0
norms_end=3142.15
FACTOR=5.011
. base_sparsifier.sh "$@"