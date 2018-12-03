#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -M alolex@vcu.edu
#PBS -N cutadapt
#PBS -j oe
#PBS -o /home/alolex/

cd /home/sequencing/work/Welm_PDX/10962R/

INFILE='/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/Welm_input01_CutAdaptAndStarAlign.list'

mkdir -p 01_cutadapt

cat $INFILE | while read s1 s2
do
	base1=`basename $s1 .fastq.gz`
	base2=`basename $s2 .fastq.gz`
	cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -u 5 -U 5 -o 01_cutadapt/$base1.trimmed.fastq.gz -p 01_cutadapt/$base2.trimmed.fastq.gz $s1 $s2
done
