#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M alolex@vcu.edu
#PBS -N SplitByOrg
#PBS -j oe
#PBS -o SplitByOrg-whim2-082317.log

cd /home/alolex/data/Harrell/WHIM2_altPipeline_8-23-17

INPUT=/home/alolex/data/Harrell/WHIM2_altPipeline_8-23-17/whim2_bam_files_to_split.list
OUTDIR=/home/alolex/data/Harrell/WHIM2_altPipeline_8-23-17

mkdir -p $OUTDIR

echo "Starting STAR alignment on "`date`

cat $INPUT | while read file
do
	echo "Starting splitting of $file on "`date` 
	prefix=$(basename $file .bam) 

	samtools view -h $file | grep ENST | samtools view -b > $OUTDIR/$prefix.human.bam
	samtools view -h $file | grep ENSMUST | samtools view -b > $OUTDIR/$prefix.mouse.bam

	echo $prefix
	echo "Finished splitting $file on "`date`

done

echo "Completed splitting on "`date`
