#!/bin/bash

rpath=/code/s6.Predict_test_data.R
exptdir=$1
agggrp=$2
numcore=$3

input=$agggrp.morefts.input_CNN_regression_modgen.RData

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir

echo =========================================================================================
echo "S6 Start: $(date)"

Rscript $rpath -i $input -t  $numcore

cp $agggrp.output_prediction_CNN_class0.txt /results/S6_${agggrp}.output_prediction_CNN_class0.txt
cp $agggrp.output_prediction_CNN_class1.txt /results/S6_${agggrp}.output_prediction_CNN_class1.txt
cp $agggrp.output_prediction_CNN_class2.txt /results/S6_${agggrp}.output_prediction_CNN_class2.txt

echo -e S6 End: $(date) "\n"
echo =========================================================================================