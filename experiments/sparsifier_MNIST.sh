#!/bin/bash
#dataset_dir=/home/ricky/Cotter-Code
#calculate_error=/home/ricky/Cotter-Code/issvm/computeMissRate.py
dataset_dir=..
calculate_error=${dataset_dir}/issvm/computeMissRate.py

METHOD=sparsifier
DATASET=MNIST
ITERATIONS=100000
TRAIN_DATA=MNIST/MNIST.binary.train.01.txt
TEST_DATA=MNIST/MNIST.binary.test.01.txt
VAL_DATA=MNIST/MNIST.binary.test.01.txt.val.
TEST_DATA=MNIST/MNIST.binary.test.01.txt.test.
INIT_DIR=$METHOD/init
MODEL_DIR=$METHOD/model
TEST_DIR=$METHOD/test
K=0.02
norms_init=0
norms_end=1800.7
FACTOR=10
. base_sparsifier.sh "$@"