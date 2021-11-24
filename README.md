# DInoPORE: Direct detection of INOsines in native RNA with nanoPORE sequencing

## OVERVIEW
DInoPORE is a computational method to detect Adenosine-to-Inosine (A-to-I) editing sites from direct-RNA sequencing data. Editing sites identified will also have their editing rate estimated. 

This GitHub repository contains the scripts to run DInoPORE. Users who wish to train their own models using DInoPORE's architecture can also find the scripts for training (see training path section below).

A sample dataset to illustrate how the pipeline runs and the typical runtime for a small sample can be found in DInoPORE's Code Ocean capsule, DOI <Insert DOI>. This dataset contains 53,283 direct-RNA-seq reads from Xenopus Laevis stage 9 rep 1. R Markdown files reproducing key results from the DInoPORE manuscript can also be found in the Code Ocean capsule under /data/Rmarkdown.

Listed below are the computational environment and software used. Usage documentation can be found in the last section of this readme.

=====================================================================================

## REQUISITE SOFTWARE

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
- keras			>= 2.3.0
- multiROC		1.1.1
- optparse		1.6.6
- pacman		0.5.1
- plyr			1.8.6
- pracma		2.3.3
- scales		1.1.1
- tensorflow	>= 2.3.0
- tidyverse		1.3.1
- usefun		0.4.8
- zoo			1.8.9

DInoPORE has been tested on CentOS Linux 7 and Ubuntu 20.04.

=====================================================================================

## INSTALLATION - Create conda env

conda create -n myenv python=3.8.5 h5py=2.10.0 nanopolish=0.11.1 pillow=8.3.1 pyyaml=5.4.1 requests=2.26.0 samtools=1.9 scipy=1.7.1 

conda activate myenv
	
=====================================================================================

## USAGE

bash mainscript.sh -e <path/to/exptdir> -r <path/to/ref.fa> -n <num_threads> -g <aggregation_Group> -c <class_Reference>

	Arguments required
		-e Full path to directory containing "fast5" folder. User must have write permission for the parent directory of this directory too.

		-r Full path to reference genome FASTA

		-n Number of threads available to DInoPORE

		-g User-defined group name used to specify runs belonging to same group. Affects aggregation step when aggregating across multiple experiment runs.

		-c Class and edit rate reference containing ground truth information about a coordinates' class and edit rate
	
	Optional argument
		-d [y/n] Delete base-called fastq files? Default value is n.


E.g. bash mainscript.sh -e /data/xen_s9_r1_50k -r /data/reference/xenLae2.fa -n 15 -g xen50k -c /data/xen_s9_r1_50k/groundtruth_classification.txt

**NOTE**: mainscript.sh expects to find "fast5" directory within exptdir:
###
	path/to/exptdir
	└── fast5

=====================================================================================

## DOCUMENTATION

### Run Mainscript1.sh (Steps 1 to 3) on individual sequencing runs
	
#### (1) Basecall fast5 -> map to genome reference -> run nanopolish to extract signal
	Script(s):
	S1.Basecall_map_nanopolish.sh (input: $exptdir $ref $numcore)
	
	Output:
	${exptdir}/out_fastq_bam/${expt}.combined.fastq
	${exptdir}/out_fastq_bam/${expt}.combined.gm2.sam
	${exptdir}/out_fastq_bam/${expt}.sorted.bam
	${exptdir}/out_nanopolish/$expt.gm2.nanopolish_eventAlignOut.txt


#### (2) Process data: convert bam file to tsv and combine nanopolish into single signal for each 5-mer of a read
	Script(s):
	S2.Process_bam_nnpl.sh (input: $exptdir $ref $numcore)
	└── s2.Sam2tsv_processtsv.sh
	└── s2.Combine_raw_nnpl.sh
		└── s2.Combine_raw_nanopolish.R
	
	Output:
	${exptdir}/raw_features/$expt.tsv.txt
	${exptdir}/raw_features/$expt.gm2.nanopolish_eventAlignOut_combined.txt


#### (3) Combine bam-tsv file and combined nanopolish file to generate raw features table
	Script(s):
	S3.Generate_raw_features.sh (input: $exptdir $numcore $agggrp)
	└── S3.Generate_raw_features.sh
		└── s3.Generate_raw_feature_table.R

	Output:
	$(dirname $exptdir)/aggregate_reads/${expt}.tsv_nnpl_inAE.txt_grp${agggrp}

### Run Mainscript2.sh (Steps 4 to 6)

#### (4) Aggregate features of reads into positions
	Script(s):
	S4.Aggregate_reads_into_pos.sh (input: $exptdir $numcore $agggrp)
		└── s4.Aggregating_reads_pos.R
		
	Output:
	$(dirname $exptdir)/matrix_CNN/$agggrp.Agg.morefts.10bin.inML.txt

**NOTE**: This script aggregates all files in matrix_CNN that were processed with the same user-defined group name (-g) and generates one file per group for downstream processing.


### **From step 5 onwards, there are 2 paths: training and testing.**

### (I) Testing path: 
We provided 3 trained models for testing. Users can used our models to detect Inosine and quantify editing rate on their own data.

**NOTE**: Models 1 and 2 was trained using H9 and Xenopus Laevis data. Model 3 was trained using editing sites in H9, Xenopus Laevis, HCT116, Mus musculus (Mouse) and synthetic RNAs.

#### (5) Transform 1D into 2D data + Label data (class 0, 1 and 2 for unmodified, Inosine and SNP AG)
	Script(s):
	S5.Transform_dim.sh (input:  $exptdir $numcore $agggrp $classref)
		└── s5.Preprocess_data_matrix_inputCNN.R
		
	Output:
	$(dirname $exptdir)/matrix_CNN/$agggrp.morefts.input_CNN_regression_modgen.RData
	
	
#### (6) Predict testing data using Dinopore models + Plot ROC and PR curves for the result
	Script(s):
	S6.Predict.sh (input: $exptdir $agggrp $numcore $classref)
		└── s6.Predict_test_data.R
	
	Output:
	$(dirname $exptdir)/matrix_CNN/$agggrp.output_prediction_CNN_class0.txt
	$(dirname $exptdir)/matrix_CNN/$agggrp.output_prediction_CNN_class1.txt
	$(dirname $exptdir)/matrix_CNN/$agggrp.output_prediction_CNN_class2.txt


### (II) Training path: 

Users can train models using their own data, based on our models' architecture. Depending on the dataset size, number of epochs, batch size and GPU type, time required for training will vary.

For our data: model 3 was trained using 265,631 positions with 893,865 matrices, on GPU NVIDIA GeForce RTX3080 with batch size = 1024. Time for 1 epoch is 88 seconds and there were 900 epochs.

**NOTE 1**: User should make sure the training data set is balance between 3 classes, especially class 0 (unmod) and class 1 (Inosine or other modifications). Number of matrices for training set should be at least 50,000 for model 1 and model 2 and at least 500,000 for model 3
 
**NOTE 2**: Ground truth for training path must have 5 columns: contig, position, strand (p for positive and n for negative), edit (0,1 and 2), and rate (editing rate 0-1 for class 1 and -1 for class 0 and 2) (For example, see Dinopore/code/misc/Example_ground_truth_training.txt)

**NOTE 3**: User should make sure number of epochs is at least 500 epochs for model 1 and model 2, and at least 900 epochs for model 3
	
**NOTE 4**: Scripts for training uses GPU. If users wish to use CPU for training, please remove the following lines of code from the respective scripts: s6a.Training_classification_model_3class.R (lines 58-59), s6b.Training_classification_model_2class.R (lines 58-59), s6c.Training_regression_model.R (lines 74-75).
	
**NOTE 5**: Training scripts should be run in the same directory as the input file (usually matrix_CNN). The output file will be generated in the same directory.

#### (5) Transform 1D into 2D data + Label data (class 0, 1 and 2 for unmodified, Inosine and SNP AG) + Split into training and validation/testing data set
	Script(s):
	Rscript s5.Train_preprocess_data_matrix_inputCNN_train_val.R -t $numcore -i $input -o $output -c $classref
	where
		input=${groupName}.Agg.morefts.10bin.inML.txt
		output=${groupName}
		classref="groundtruth_training_file" (see NOTE 2 above)
		
	Output:
	${groupName}.validation_matrix.rds (used as $vali in step 6a-c)
	${groupName}.training_matrix.rds (used as $train in step 6a-c)
	
#### (6a) Train 3-class classification model
	Script(s):
	Rscript s6a.Training_classification_model_3class.R -v $vali -t $train -o $groupName -e $epoch -b $batch -s $seed
	where
		vali=${groupName}.validation_matrix.rds
		train=${groupName}.training_matrix.rds
		
	Output:
	${groupName}.best_pos5_mix_3class_resnet.h5
	
#### (6b) Train 2-class classification model
	Script(s):
	Rscript s6b.Training_classification_model_2class.R  -v $vali -t $train -o $groupName -e $epoch -b $batch -s $seed
		
	Output:
	${groupName}.best_pos5_mix_3c_1vs1_resnet.h5
	
#### (6c) Train regression model
	Script(s):
	Rscript s6c.Training_regression_model.R -v $vali -t $train -o $groupName -e $epoch -b $batch -s $seed
		
	Output:
	${groupName}.best_regression_morefts_16384_1024_asin06.h5
	
#### (7) Predict testing data using trained models in step 6a, 6b and 6c
	Script(s):
	Rscript s7.Predict_test_data_using_trained_models.R -i $input -t  $numcore -m $model1 -M $model2 -r $model3
	where
		input=${groupName}.validation_matrix.rds
		model1=${groupName}.best_pos5_mix_3class_resnet.h5
		model2=${groupName}.best_pos5_mix_3c_1vs1_resnet.h5
		model3=${groupName}.best_regression_morefts_16384_1024_asin06.h5

	
	Output:
	${groupName}.output_prediction_CNN_class0.txt
	${groupName}.output_prediction_CNN_class1.txt
	${groupName}.output_prediction_CNN_class2.txt
	
