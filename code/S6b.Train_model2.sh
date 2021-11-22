#!/bin/bash

codedir=$(cd `dirname $0` && pwd)

rpath=${codedir}/s6b.Training_classification_model_2class.R
inputdir=$1
groupName=$2
epoch=$3
batch=$4
seed=$5

vali=$groupName.validation_matrix.rds
train=$groupName.training_matrix.rds

cd $inputdir

echo =========================================================================================
echo "S6b Start: $(date)"

Rscript $rpath -v $vali -t $train -o $groupName -e $epoch -b $batch -s $seed

echo -e "S6b End: $(date) \n"
echo =========================================================================================