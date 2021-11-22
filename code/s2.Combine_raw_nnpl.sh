#!/bin/bash

npdir=$1
expt=$2
exptdir=$3
numcore=$4

codedir=$(cd `dirname $0` && pwd)

cd $npdir

rpath=${codedir}/s2.Combine_raw_nanopolish.R
file=${exptdir}/out_nanopolish/${expt}.gm2.nanopolish_eventAlignOut.txt
head -n 1 $file > raw_nanopolish.header
noline=$(wc -l $file | cut -d " " -f 1)
out=$(echo $file | sed -e "s/.txt/_combined.txt/g")
tmp=tmp.nnpl

echo $file
echo $noline
echo $out

num=2
chunk=10000000
num1=$(expr $num + $chunk)
num2=$(expr $num1 + 1)
filecount=1

while [ $num -le $noline ]
do

	sed -n "$num,${num1}p; ${num2}q" $file > $tmp
	chk=$(wc -l $tmp | cut -f 1 -d " ")
	echo Processing reads $num to $num1 in file $filecount

	if [ $chk -gt 0 ]; then
		Rscript $rpath -f $tmp -t $numcore -o ${out}.part${filecount} -s $noline -n $num -c $chunk
		
		if test -f tmp.eli; then
			eli=$(cat tmp.eli)
			rm tmp.eli
		else
			eli=0
		fi
				
	fi
	
	num=$(( $num1 - $eli + 1 ))
	num1=$(expr $num + $chunk)
	num2=$(expr $num1 + 1)
	filecount=$(expr $filecount + 1)

done

cat ${codedir}/misc/nnpl.header > $out
cat ${out}.part* | grep -v contig >> $out

rm tmp.nnpl ${out}.part*

ftdir=${exptdir}/raw_features
mv $out -t $ftdir
