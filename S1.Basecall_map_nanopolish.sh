#!/bin/bash
 
flowcell=FLO-MIN106
kit=SQK-RNA002

#Print run details and define variable names
echo =========================================================================================
echo "S1 Start: $(date)"
echo Flowcell used: $flowcell
echo Sequencing kit used: $kit

#Name variables
exptdir=$1 #directory of sample folder. For example: ../raw_data/<name expt>
ref=$2
numcore=$3
delfastq=$4

expt=$(basename $exptdir)
fast5dir=${exptdir}/fast5
fqdir=${exptdir}/out_fastq_bam
npdir=${exptdir}/out_nanopolish
fqCombined=${expt}.combined.fastq
fq=$fqdir/$fqCombined
sam=${expt}.combined.gm2.sam
sortedBam=$expt.sorted.bam
npsummary=$expt.$(basename $ref).gm2.nanopolish.sum
np_out=$npdir/$expt.gm2.nanopolish_eventAlignOut.txt

cd $exptdir
mkdir -p $fqdir $npdir

#Step 1 - guppy; basecall fastq from fast5/
#guppy_basecaller -i $fast5dir -s $fqdir --flowcell $flowcell --kit $kit --cpu_threads_per_caller 12 -r

#Step 2 - Concatenate all fastq files into one main fastq file "*.combined.fastq".
cd $fqdir
#rm $fqCombined
#cat *.fastq > $fqCombined
#echo All fastq have been combined.

#if [ "$delfastq" = "y" ]; then
#  find . -maxdepth 1 -type f -name '*.fastq' -a ! -name '*.combined.fastq' -delete
#elif [ "$delfastq" = "n" ]; then
#  echo "Base-called FASTQ files preserved. Proceeding to map reads."
#fi																					

#Step 3 - Map to reference genome. Generates sam file.
sed -i 's/U/T/g' $fqCombined
$graphmap2 align -x rnaseq -t $numcore -r $ref -d $fqCombined -o $sam
echo "graphmap2 has finished aligning samples to the reference with 'graphmap2_nopost align -x rnaseq -t $numcore -r $ref -d $fqCombined -o $sam'."

#Step 4 - Index bam and generate statistics - coverage, error rate, number of reads.
samtools view -@ $numcore -S $sam -b | samtools sort -o $sortedBam - ; samtools index $sortedBam

#Step 5 - Back in the main folder, index fast5 to fastq with nanopolish index
nanopolish index -d $fast5dir/ $fq

#Step 6 - Run nanopolish for fastq files to compute an improved consensus sequence from the draft genome assembly produced by minimap2
nanopolish eventalign --reads $fq --bam $sortedBam --genome $ref --summary $npsummary --print-read-names --threads $numcore --scale-events > $np_out

echo -e S1 End: $(date) "\n"
echo =========================================================================================