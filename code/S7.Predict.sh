#!/bin/bash

rpath=/code/s7.Predict_test_data_using_trained_models.R
exptdir=$1
agggrp=$2
numcore=$3
input=$4
model1=$5
model2=$6
model3=$7

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir

echo =========================================================================================
echo "S7 Start: $(date)"

Rscript $rpath -i $input -t  $numcore -m3 $model1 -m2 $model2 -r $model3

cp $agggrp.output_prediction_CNN_class0.txt /results/S7_${agggrp}.output_prediction_CNN_class0.txt
cp $agggrp.output_prediction_CNN_class1.txt /results/S7_${agggrp}.output_prediction_CNN_class1.txt
cp $agggrp.output_prediction_CNN_class2.txt /results/S7_${agggrp}.output_prediction_CNN_class2.txt

echo -e S7 End: $(date) "\n"
echo =========================================================================================