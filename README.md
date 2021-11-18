# Dinopore
The scripts used to generate the results found in our manuscript can be found in this capsule. Also included in the capsule is a sample dataset containing 53,283 direct-RNA-seq reads from Xenopus Laevis to illustrate how the pipeline runs and the typical runtime for a small sample (Just hit "Reproducible Run" on the right!). 

Rmarkdown files reproducing results from the manuscript can also be found in /data/Rmarkdown.

Listed below are the computational environment and software used. Usage documentation can be found in the last section of this readme.

Capsule DOI: <Insert DOI>

Environment build time is approximately 15 minutes.

=====================================================================================

REQUISITE SOFTWARE

Guppy_basecaller 3.2.4 (Please modify S1.Basecall_map_nanopolish.sh to point to it)

Third party software included
- Graphmap2		 			0.6.3 
- Sam2tsv (part of Jvarkit)	34d8e7f7 
- Picard					2.21.6

conda env
- python		3.8.5
- h5py			2.10.0
- nanopolish	0.11.1
- pillow		8.3.1
- pyyaml		5.4.1
- requests		2.26.0
- samtools		1.9
- scipy			1.7.1
- tensorflow	2.4.1

R Packages (R version 4.1)
- Matrix		1.3.4
- R.utils		2.11.0
- Rcpp			1.0.7
- abind			1.4.5
- caret			6.0.89
- data.table	1.14.2
- doParallel	1.0.16
- ff			4.0.4
- foreach		1.5.1
- keras			2.6.1
- multiROC		1.1.1
- optparse		1.6.6
- pacman		0.5.1
- plyr			1.8.6
- pracma		2.3.3
- scales		1.1.1
- tensorflow	2.6.0
- tidyverse		1.3.1
- usefun		0.4.8
- zoo			1.8.9

DInoPORE has been tested on CentOS Linux 7 and Ubuntu 20.04.

=====================================================================================

INSTALLATION

conda create -n myenv python=3.8.5 h5py=2.10.0 nanopolish=0.11.1 pillow=8.3.1 pyyaml=5.4.1 requests=2.26.0 samtools=1.9 scipy=1.7.1 tensorflow=2.4.1

conda activate myenv

=====================================================================================

USAGE

bash mainscript.sh -e <path/to/exptdir> -r <path/to/ref.fa> -n <num_threads> -g <aggregationGroup> -c <classReference>

	Arguments required
		-e Full path to directory containing "fast5" folder. User must have write 								permission for the parent directory of this directory too.

		-r Full path to reference genome FASTA

		-n Number of threads available to DInoPORE

		-g User-defined group name used to specify runs belonging to same group. Affects 						aggregation step when aggregating across multiple experiment runs.

		-c Class and edit rate reference containing ground truth information about a 							coordinates' class and edit rate
	
	Optional argument
		-d [y/n] Delete base-called fastq files? Default value is n.


E.g. bash mainscript.sh -e /data/xen_s9_r1_50k -r /data/reference/xenLae2.fa -n 15 -g xen50k -c /data/xen_s9_r1_50k/groundtruth_classification.txt

Note: mainscript.sh expects to find "fast5" directory within exptdir
i.e. 
	path/to/exptdir
	└── fast5

=====================================================================================

DOCUMENTATION
Steps:

(1) Basecall fast5 -> map to genome reference -> run nanopolish to extract signal
	Script:
	S1.Basecall_map_nanopolish.sh (input: $exptdir $ref $numcore)
	
	Output:
	${exptdir}/out_fastq_bam/${expt}.combined.fastq
	${exptdir}/out_fastq_bam/${expt}.combined.gm2.sam
	${exptdir}/out_fastq_bam/${expt}.sorted.bam
	${exptdir}/out_nanopolish/$expt.gm2.nanopolish_eventAlignOut.txt


(2) Process data: convert bam file to tsv and combine nanopolish into single signal for each 5-mer of a read
	Scripts:
	S2.Process_bam_nnpl.sh (input: $exptdir $ref $numcore)
	└── s2.Sam2tsv_processtsv.sh
	└── s2.Combine_raw_nnpl.sh
		└── s2.Combine_raw_nanopolish.R
	
	Output:
	${exptdir}/raw_features/$expt.tsv.txt
	${exptdir}/raw_features/$expt.gm2.nanopolish_eventAlignOut_combined.txt


(3) Combine bam-tsv file and combined nanopolish file to generate raw features table
	Scripts:
	S3.Generate_raw_features.sh (input: $exptdir $numcore $agggrp)
	└── S3.Generate_raw_features.sh
		└── s3.Generate_raw_feature_table.R

	Output:
	$(dirname $exptdir)/aggregate_reads/${expt}.tsv_nnpl_inAE.txt_grp${agggrp}


(4) Aggregate features of reads into positions
	Scripts:
	S4.Aggregate_reads_into_pos.sh (input: $exptdir $numcore $agggrp)
		└── s4.Aggregating_reads_pos.R
		
	Output:
	$(dirname $exptdir)/matrix_CNN/$agggrp.Agg.morefts.10bin.inML.txt

	Note: This script can be run on any number of runs and will determine how reads across runs (if applicable) are aggregated before passing into the CNN model.

(5) Transform 1D into 2D data + Label data (class 0, 1 and 2 for unmodified, Inosine and SNP AG)
	Scripts:
	S5.Transform_dim.sh (input:  $exptdir $numcore $agggrp $classref)
		└── s5.Preprocess_data_matrix_inputCNN.R
		
	Output:
	$(dirname $exptdir)/matrix_CNN/$agggrp.morefts.input_CNN_regression_modgen.RData
	
	
(6) Predict testing data using Dinopore models + Plot ROC and PR curves for the result
	S6.Predict.sh (input: $exptdir $agggrp $numcore $classref)
		└── s6.Predict_test_data.R
	
	Output:
	$(dirname $exptdir)/matrix_CNN/$agggrp.output_prediction_CNN_class0.txt
	$(dirname $exptdir)/matrix_CNN/$agggrp.output_prediction_CNN_class1.txt
	$(dirname $exptdir)/matrix_CNN/$agggrp.output_prediction_CNN_class2.txt
