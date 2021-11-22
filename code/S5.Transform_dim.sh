#!/bin/bash

exptdir=$1
numcore=$2
agggrp=$3
classref=$4

codedir=$(cd `dirname $0` && pwd)

input=$agggrp.Agg.morefts.10bin.inML.txt
output=$agggrp.morefts.input_CNN_regression_modgen.RData

rpath=${codedir}/s5.Preprocess_data_matrix_inputCNN.R

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir

echo =========================================================================================
echo "S5 Start: $(date)"

Rscript $rpath -t $numcore -i $input -o $output -c $classref 

echo -e S5 End: $(date) "\n"
echo =========================================================================================
