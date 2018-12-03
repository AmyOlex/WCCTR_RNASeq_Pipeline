#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M alolex@vcu.edu
#PBS -N alignStats
#PBS -j oe
#PBS -o alignStats-official-081717.log

cd /home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow

INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/input03_getAlignStats.list
OUTFILE=/home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/output03_bamStats_081717.txt
UMAP1=".Unmapped.out.mate1"
UMAP2=".Unmapped.out.mate2"


touch $OUTFILE

echo "Starting file counting on "`date`
echo -e "Sample,Total Counts,Human Aligned,Mouse Aligned,Unmapped\n" >> $OUTFILE

cat $INPUT | while read file
do
	echo "Starting processing of $file on "`date` 
	dir=`dirname $file`
	base=`basename $file .Aligned.out.bam`
	umapfile1=$dir/$base$UMAP1
	umapfile2=$dir/$base$UMAP2

	
	subtotal=`samtools view $file | wc -l`
	human=`samtools view $file | grep HUMAN | wc -l`
	mouse=`samtools view $file | grep MOUSE | wc -l`
	u1=`cat $umapfile1 | grep @ | wc -l`
	u2=`cat $umapfile2 | grep @ | wc -l`
	total=$(($subtotal+$u1+$u2))
	utotal=$(($u1+$u2))

	echo -e "$base,$total,$human,$mouse,$utotal\n" >> $OUTFILE
	echo "Finished processing of $file on "`date`

done

echo "Completed bam processing on "`date`
