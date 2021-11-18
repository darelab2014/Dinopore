#!/bin/bash

#Note: remember to unzip fasta
#This script is used to convert bam file into tsv file, then process it (remove S, H and N, only keep M (match or mismatch), D (deletion) and I (insertion). Then replace NA in called base, format positive and negative annotation)

fqdir=$1
expt=$2
exptdir=$3
ref=$4
numcore=$5

cd $fqdir

sortedBam=$expt.sorted.bam
final=$expt.tsv.txt
ftdir=${exptdir}/raw_features

mkdir -p $ftdir

#Process tsv file to remove S, H and N, only keep M (match or mismatch), D (deletion) and I (insertion). Then replace NA in called base, format positive and negative annotation
samtools view -@ $numcore -h -F 4 $sortedBam | $SAM2TSV -r $ref | awk 'BEGIN{FS=OFS="\t"} ($9 != "S") && ($9 != "H") && ($9 != "N")' - | awk 'BEGIN{FS=OFS="\t"} ($7=="."){$7="-99";} ($4=="."){$4="-99"} ($5=="."){$5="na"} ($8=="."){$8="na"} ($9=="D"){$6=" "} ($2==16){$2="n"} ($2==0){$2="p"} 1' > $final

mv $final -t $ftdir
