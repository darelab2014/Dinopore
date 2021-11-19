#!/bin/bash

rpath=/code/s6b.Training_classification_model_2class.R
exptdir=$1
agggrp=$2
numcore=$3
vali=$agggrp.validation_matrix.rds
train=$agggrp.training_matrix.rds
epoch=$4
batch=$5
seed=$6

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir
echo =========================================================================================
echo "S6b Start: $(date)"

Rscript $rpath -v $vali -t $train -o $agggrp -e $epoch -b $batch -s $seed

echo -e S6b End: $(date) "\n"
echo =========================================================================================