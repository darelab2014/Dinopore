#!/bin/bash

rpath=/code/s5.Preprocess_data_matrix_inputCNN.R
exptdir=$1
numcore=$2
agggrp=$3
classref=$4
input=$agggrp.Agg.morefts.10bin.inML.txt
output=$agggrp.morefts.input_CNN_regression_modgen.RData

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir

echo =========================================================================================
echo "S5 Start: $(date)"

Rscript $rpath -t $numcore -i $input -o $output -c $classref 

#cp $output /results/S5_$output

echo -e S5 End: $(date) "\n"
echo =========================================================================================