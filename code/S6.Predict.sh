#!/bin/bash

exptdir=$1
agggrp=$2
numcore=$3

codedir=$(cd `dirname $0` && pwd)

input=$agggrp.morefts.input_CNN_regression_modgen.RData

rpath=${codedir}/s6.Predict_test_data.R

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir

echo =========================================================================================
echo "S6 Start: $(date)"

Rscript $rpath -i $input -t  $numcore

echo -e S6 End: $(date) "\n"
echo =========================================================================================