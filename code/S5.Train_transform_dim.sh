#!/bin/bash

codedir=$(cd `dirname $0` && pwd)

rpath=${codedir}/s5.Train_preprocess_data_matrix_inputCNN_train_val.R
inputdir=$1
numcore=$2
groupName=$3
classref=$4

input=$groupName.Agg.morefts.10bin.inML.txt
output=$groupName

cd $inputdir

echo =========================================================================================
echo "S5 Start: $(date)"

Rscript $rpath -t $numcore -i $input -o $output -c $classref 

echo -e "S5 End: $(date) \n"
echo =========================================================================================