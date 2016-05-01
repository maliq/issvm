#!/bin/bash
#dataset_dir=/home/ricky/Cotter-Code
#calculate_error=/home/ricky/Cotter-Code/issvm/computeMissRate.py
dataset_dir=..
calculate_error=${dataset_dir}/issvm/computeMissRate.py

METHOD=sparsifier
DATASET=Adult
ITERATIONS=100000
TRAIN_DATA=Adult/a8a.original.train.01.txt
VAL_DATA=Adult/a8a.original.test.01.txt.val.
TEST_DATA=Adult/a8a.original.test.01.txt.test.
INIT_DIR=$METHOD/init
MODEL_DIR=$METHOD/model
TEST_DIR=$METHOD/test
K=0.1
norms_init=0
norms_end=1061.1
#factor to split norm^2
FACTOR=10
. base_sparsifier.sh "$@"
