#!/bin/bash
#dataset_dir=/home/ricky/Cotter-Code
#calculate_error=/home/ricky/Cotter-Code/issvm/computeMissRate.py
dataset_dir=..
calculate_error=${dataset_dir}/issvm/computeMissRate.py

METHOD=sparsifier
DATASET=IJCNN
ITERATIONS=100000
TRAIN_DATA=IJCNN/ijcnn1.original.train.01.txt
TEST_DATA=IJCNN/ijcnn1.original.test.01.txt
VAL_DATA=IJCNN/ijcnn1.original.test.01.txt.val.
TEST_DATA=IJCNN/ijcnn1.original.test.01.txt.test.
INIT_DIR=$METHOD/init
MODEL_DIR=$METHOD/model
TEST_DIR=$METHOD/test
K=1.0
declare -a norms=(38.72 71.47 131.61 230.26 421.38 880.83 2074.61 5043.84 10305.4 12044.6)
norms_init=38.72
norms_end=0
. base_sparsifier.sh ${1} ${2}