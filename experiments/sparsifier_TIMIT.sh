#!/bin/bash
#dataset_dir=/home/ricky/Cotter-Code
#calculate_error=/home/ricky/Cotter-Code/issvm/computeMissRate.py
dataset_dir=..
calculate_error=${dataset_dir}/issvm/computeMissRate.py

METHOD=sparsifier
DATASET=TIMIT
ITERATIONS=100000
TRAIN_DATA=TIMIT/TIMIT.3.binary.train.01.txt
TEST_DATA=TIMIT/TIMIT.3.binary.test.01.txt
VAL_DATA=TIMIT/TIMIT.3.binary.test.01.txt.val.
TEST_DATA=TIMIT/TIMIT.3.binary.test.01.txt.test.
INIT_DIR=$METHOD/init
MODEL_DIR=$METHOD/model
TEST_DIR=$METHOD/test
K=0.025
#declare -a norms=(24.13 8050.47)
norms_init=0
norms_end=8050.47
FACTOR=10
. base_sparsifier.sh "$@"