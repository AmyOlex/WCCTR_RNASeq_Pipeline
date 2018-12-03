#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M alolex@vcu.edu
#PBS -N insertSize
#PBS -j oe
#PBS -o insertSize.log

cd /home/sequencing/data/WorkData/Turner-Harrell_SCS_AllBatches

#INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/Turner_SCS_input06_getInsertSize.list
INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/Turner_SCS_input06_getInsertSizeTranscriptome.list
OUTDIR=/home/sequencing/data/WorkData/Turner-Harrell_SCS_AllBatches/06_getInsertSizeTranscriptome

mkdir -p $OUTDIR

cat $INPUT | while read file
do
    echo "Starting $file"
    base=`basename $file Aligned.toTranscriptome.out.bam`
    samtools view -s .10 $file | cut -f 9 | sort -n > $OUTDIR/$base.insertSizes.txt
    echo "Done $file"
done
