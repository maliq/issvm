#!/usr/bin/env bash
declare -a DS_LIST=(Adult IJCNN MNIST TIMIT Web)
for DS in "${DS_LIST[@]}"
do
    ln -s ../issvm/experiments/sparsifier_${DS}.sh sparsifier_${DS}.sh
done
ln -s ../issvm/experiments/base_sparsifier.sh base_sparsifier.sh