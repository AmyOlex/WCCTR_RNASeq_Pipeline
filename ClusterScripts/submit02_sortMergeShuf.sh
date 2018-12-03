#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M alolex@vcu.edu
#PBS -N SMS
#PBS -j oe
#PBS -o sort-merge-shuf-official-081717.log

cd /home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow

INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/input02_sort-merge-shuf.list
OUTDIR=/home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/01_star-align
TMPDIR=/home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/02_sort-merge-tmp

mkdir -p $TMPDIR
mkdir -p $OUTDIR

echo "Starting sort-merge-shuffle on "`date`

cat $INPUT | while read file1 file2
do
	echo "Starting sort of $file1 and $file2 on "`date` 
	prefix1=$(basename $file1 001.Aligned.toTranscriptome.out.bam) 
	prefix2=$(basename $file2 001.Aligned.toTranscriptome.out.bam)

	samtools sort -@ 12 -o $TMPDIR/$prefix1.sort.bam $file1
	samtools sort -@ 12 -o $TMPDIR/$prefix2.sort.bam $file2
	
	samtools merge -r $TMPDIR/$prefix1-$prefix2.merged.bam $TMPDIR/$prefix1.sort.bam $TMPDIR/$prefix2.sort.bam 

	samtools collate $TMPDIR/$prefix1-$prefix2.merged.bam $OUTDIR/$prefix1-$prefix2.merged.shuf.Aligned.toTranscriptome.out

	echo "Finished merging $prefix1 and $prefix2 on "`date`

done

echo "Completed sort-merge-shuffle on "`date`
