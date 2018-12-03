#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M alolex@vcu.edu
#PBS -N fastqc
#PBS -j oe


cd /home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow

INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/input00_fastqc.list

mkdir -p 00_fastqc

OUTDIR=/home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/00_fastqc

cat $INPUT | while read file
do 
	echo $file;
	/home/alolex/bin/FastQC/fastqc -t 4 -o $OUTDIR --noextract $file;
done
