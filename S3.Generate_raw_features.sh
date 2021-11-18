#!/bin/bash

exptdir=$1
numcore=$2
agggrp=$3

rpath=/code/s3.Generate_raw_feature_table.R
expt=$(basename ${exptdir})
ftdir=${exptdir}/raw_features
cd $ftdir

echo =========================================================================================
echo "S3 Start: $(date)"

#Name variables
nnpl=$expt.gm2.nanopolish_eventAlignOut_combined.txt
tsv=${expt}.tsv.txt
outname=${expt}.tsv_nnpl_inAE.txt
intsv=${expt}.input.filteredSHN.tsv.txt
innnpl=${expt}.input.nanopolish_eventAlignOut_combined.txt
readtsv=${expt}.read.tsv
readnnpl=${expt}.read.nnpl
interread=${expt}.intersect.read

#Step 1 - Extract common reads between tsv and nanopolish file
awk 'NR!=1 {print $1}' $tsv | LC_ALL=C uniq | LC_ALL=C sort -S 80% --parallel=$numcore  > $readtsv
awk '{print $2}' $nnpl | LC_ALL=C uniq | LC_ALL=C sort -S 80% --parallel=$numcore  > $readnnpl
comm -12 $readtsv $readnnpl > $interread
	
num=1
chunk=100000
num1=$(expr $num + $chunk - 1)
num2=$(expr $num1 + 1)
noline=$(wc -l $interread | cut -d " " -f 1)
filecount=1

while [ $num -le $noline ]
do
	echo Processing reads $num to $num1
	
	sed -n "$num,${num1}p; ${num2}q" $interread > readstmp
	chk=$(wc -l readstmp | cut -f 1 -d " ")
	
	if [ $chk -gt 0 ]; then
		cat /code/misc/nnpl.header > $innnpl
		cat /code/misc/tsv.header > $intsv

		LC_ALL=C grep -h -F -w -f readstmp $tsv >> $intsv
		LC_ALL=C grep -h -F -w -f readstmp $nnpl >> $innnpl

		#Step 2 - Generate raw features files
		Rscript $rpath -n $innnpl -t $intsv -o ${outname}.part$filecount
												   

		num=$num2
		num1=$(expr $num + $chunk - 1)
		num2=$(expr $num1 + 1)
		filecount=$(expr $filecount + 1)
		
		rm readstmp
	fi	
done

cat /code/misc/inae.header > $outname
cat ${outname}.part*.positive | LC_ALL=C sort -k1,1 -k3,3n | grep -v contig >> $outname
cat ${outname}.part*.negative | LC_ALL=C sort -k1,1 -k3,3n | grep -v contig >> $outname

echo ${outname} generated on $(date)

rm $readtsv $readnnpl $interread $intsv $innnpl ${outname}.part*.positive ${outname}.part*.negative

aggdir=$(dirname $exptdir)/aggregate_reads
mkdir $aggdir
#cp $outname /results/S3_${outname}_grp${agggrp}
mv $outname ${outname}_grp${agggrp}
mv ${outname}_grp${agggrp} -t $aggdir

echo -e S3 End: $(date) "\n"
echo =========================================================================================