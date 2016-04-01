#!/bin/bash
dataset_dir=/home/ricky/Cotter-Code
calculate_error=/home/ricky/Cotter-Code/issvm/computeMissRate.py

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
declare -a norms=(24.13 8050.47)
norms_init=2.15
norms_end=1800.7
. base_sparsifier.sh ${1} ${2}