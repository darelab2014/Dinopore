#!/bin/bash

#This script with s4.Aggregating_reads_pos.R are used to aggregate features of all reads of the same positions
#Output: Agg.10bin.inML.txt file: aggreagated file. Output file will be moved to next step folder (work/matrix_CNN)


codedir=$(cd `dirname $0` && pwd)

exptdir=$1
numcore=$2
agggrp=$3
outname=$agggrp.Agg.morefts.10bin.inML.txt
rpath=${codedir}/s4.Aggregating_reads_pos.R

cd $(dirname $exptdir)/aggregate_reads
cnndir=$(dirname $exptdir)/matrix_CNN
mkdir -p $cnndir

echo =========================================================================================
echo "S4 Start: $(date)"

Rscript $rpath -t $numcore -o $outname -r $agggrp

mv $outname -t ${cnndir}

echo -e S4 End: $(date) "\n"
echo =========================================================================================
