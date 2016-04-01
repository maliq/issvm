#!/bin/bash
dataset_dir=/home/ricky/Cotter-Code
calculate_error=/home/ricky/Cotter-Code/issvm/computeMissRate.py

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
#declare -a norms=(15.01 25.85 45.99 85.51 158.86 280.88 457.40 685.79 934.99 1061.1)
#declare -a norms=(15.01 25.85 45.99 85.51 158.86 219.87 280.88 457.40 685.79 934.99 1061.1)
#NORM=8.95435e-08
norms_init=15.01
norms_end=1061.1
. base_sparsifier.sh ${1} ${2}