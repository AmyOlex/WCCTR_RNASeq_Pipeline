#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M alolex@vcu.edu
#PBS -N countClusters
#PBS -j oe

cd /home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/

INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/Harrell_All_input00a_fastqClusterCounts2.list

OUTFILE=/home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/output_00a_fastqClusterCountsAll2.txt

touch $OUTFILE

cat $INPUT | while read file
do 
	counts=`expr $(zcat $file | wc -l) / 4` >> $OUTFILE
	echo -e $file "\t" $counts >> $OUTFILE
done

