#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M alolex@vcu.edu
#PBS -N insertSize
#PBS -j oe
#PBS -o insertSize.log

cd /home/sequencing/data/WorkData/Turner-Harrell_SCS_AllBatches

#INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/Turner_SCS_input07_averageInsertSize.list
INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/Turner_SCS_input07_averageInsertSizeTranscriptome.list
OUTDIR=/home/sequencing/data/WorkData/Turner-Harrell_SCS_AllBatches/06_getInsertSizeTranscriptome

mkdir -p $OUTDIR

cat $INPUT | while read file
do
    echo "Starting $file"
    base=`basename $file .insertSizes.txt`
    lines=$(cat $file | wc -l)
    (( e=$lines / 2 ))
    tail -n $e $file | grep [^0] | awk '{ total += $1; count++ } END { print total/count }' > $OUTDIR/$base.average.txt
    
    echo "Done $file"
done
