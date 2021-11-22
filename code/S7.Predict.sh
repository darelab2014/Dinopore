#!/bin/bash

codedir=$(cd `dirname $0` && pwd)

rpath=${codedir}/s7.Predict_test_data_using_trained_models.R
inputdir=$1
numcore=$2
groupName=$3

input=$groupName.validation_matrix.rds

model1=$groupName.best_pos5_mix_3class_resnet.h5
model2=$groupName.best_pos5_mix_3c_1vs1_resnet.h5
model3=$groupName.best_regression_morefts_16384_1024_asin06.h5

cd $inputdir

echo =========================================================================================
echo "S7 Start: $(date)"

Rscript $rpath -i $input -t  $numcore -m $model1 -M $model2 -r $model3

echo -e "S7 End: $(date) \n"
echo =========================================================================================