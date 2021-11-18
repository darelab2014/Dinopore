#!/bin/bash

exptdir=$1
ref=$2
numcore=$3
expt=$(basename ${exptdir})
fqdir=${exptdir}/out_fastq_bam
npdir=${exptdir}/out_nanopolish
dict=$(echo $ref | sed 's/.fasta//g' | sed 's/.fa//g').dict

echo =========================================================================================
echo "S2 Start: $(date)"

#Part 1 - sam2tsv
#Step 1 - Create Sequence Dictionary
if [ ! -f "$dict" ]; then
	echo Sequence dictionary not found. Creating sequence dictionary.
	$PICARD CreateSequenceDictionary R=$ref O=$dict
else
	echo Sequence dictionary found.
fi

#Step 2 - Run sam2tsv and process tsv file
sh ./s2.Sam2tsv_processtsv.sh $fqdir $expt $exptdir $ref $numcore

#Part 2 - process raw nanopolish file
sh ./s2.Combine_raw_nnpl.sh $npdir $expt $exptdir $numcore

echo -e S2 End: $(date) "\n"
echo =========================================================================================