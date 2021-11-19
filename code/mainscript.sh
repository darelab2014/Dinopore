#!/usr/bin/env bash
set -ex

# This is the master script for the capsule. When you click "Reproducible Run", the code in this file will execute.
								

export SAM2TSV='java -jar /code/misc/sam2tsv.jar'
export PICARD='java -jar /code/misc/picard.jar'
export graphmap2='/code/misc/graphmap2'

chmod 744 $graphmap2
delfastq=n

usage() {                                 # Function: Print a help message.
    echo ""
    echo "Usage: $0 [ -e PATH/TO/EXPTDIR ] [ -r PATH/TO/REFERENCE.FASTA ] [ -n num_threads ] [ -g aggregateGroup ] [ -c classReference ]" 1>&2
    echo ""
    echo "Optional:"
    echo "  -d  [y/n] Delete base-called fastq files? Default value is n."
    echo ""
}
exit_abnormal() {                         # Function: Exit with error.
  usage
  exit 1
}
while getopts ":e:r:n:g:c:d:" options; do         # Loop: Get the next option;
                                          # use silent error checking;
                                          # options n and t take arguments.
  case "${options}" in                    #
    e)                                    # If the option is n,
      exptdir=${OPTARG}                   # set $exptdir to specified value.
      ;;
    r)                                    # If the option is t,
      ref=${OPTARG}                       # Set $ref to specified value.
      ;;
    n)
      numcore=${OPTARG}
      re_isanum='^[0-9]+$'                # Regex: match whole numbers only
      if ! [[ $numcore =~ $re_isanum ]] ; then   # if $numcore not whole:
        echo "Error: numcore must be a positive integer."
        exit_abnormal
        exit 1
      elif [ $numcore -eq "0" ]; then     # If it's zero:
        echo "Error: numcore must be a positive integer"
        exit_abnormal                     # Exit abnormally.
      fi
      ;;
    g)                                    # If unknown (any other) option:
	  agggrp=${OPTARG}					  # Set $agggrp to specified value.
      ;;
    c)                                   # If unknown (any other) option:
	  classref=${OPTARG}				  # Set $classref to specified value.
      ;;
    d)                                   # If unknown (any other) option:
	  delfastq=${OPTARG}				  # Set $delfastq to specified value.
      ;;
	  :)                                    # If expected argument omitted:
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal                       # Exit abnormally.
      ;;
    *)                                    # If unknown (any other) option:

      exit_abnormal                       # Exit abnormally.
      ;;
  esac
done

if [ "$exptdir" = "" ]; then
        echo "Error: please provide path to experiment directory."
        exit_abnormal
fi

if [ "$ref" = "" ]; then
        echo "Error: please provide path to FASTA for mapping."
        exit_abnormal
fi

if [ "$numcore" = "" ]; then
        echo "Error: please indicate number of cores to assign to DInoPORE."
        exit_abnormal
fi

if [ "$agggrp" = "" ]; then
        echo "Error: please specify a group name for aggregation. Please use only letters and/or numbers. Avoid using symbols."
        exit_abnormal
fi

if [ "$classref" = "" ]; then
        echo "Error: please provide class and edit rate ground truth."
        exit_abnormal
        exit 1
fi

if [ "$delfastq" = "" ]; then
  delfastq=n
elif [ "$delfastq" = "n" ]; then
  echo "FASTQ files in out_fastq_bam folder will be preserved"
elif [ "$delfastq" = "N" ]; then
  delfastq=n
  echo "FASTQ files in out_fastq_bam folder will be preserved"
elif [ "$delfastq" = "y" ]; then
  echo "FASTQ files in out_fastq_bam folder will be deleted"
elif [ "$delfastq" = "Y" ]; then
  delfastq=y
  echo "FASTQ files in out_fastq_bam folder will be deleted"
else
  echo "Error: [y/N] please specify whether FASTQ files in path/to/exptdir/out_fastq_bam should be deleted. WARNING: ALL FILES ENDING WITH *.FASTQ WILL BE DELETED. IF UNSURE, INPUT "n"."
  exit_abnormal
  exit 1
fi


#Step 1 - Basecall fast5 -> map to genome reference -> run nanopolish to extract signal
sh ./S1.Basecall_map_nanopolish.sh $exptdir $ref $numcore $delfastq

#Step 2 - Process data: convert bam file to tsv and combine nanopolish into single signal for each 5-mer of a read
sh ./S2.Process_bam_nnpl.sh $exptdir $ref $numcore

#Step 3 - Combine bam-tsv file and combined nanopolish file to generate raw features table
sh ./S3.Generate_raw_features.sh $exptdir $numcore $agggrp

#Step 4 - Aggregate features of reads into positions
sh ./S4.Aggregate_reads_into_pos.sh $exptdir $numcore $agggrp

#Step 5 - Transform 1D into 2D data + Label data (class 0, 1 and 2 for unmodified, Inosine and SNP AG)
sh ./S5.Transform_dim.sh $exptdir $numcore $agggrp $classref

#Step 6 - Predict testing data using Dinopore models + Plot ROC and PR curves for the result
sh ./S6.Predict.sh $exptdir $agggrp $numcore
