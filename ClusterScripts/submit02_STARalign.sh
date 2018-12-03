#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=2:ppn=24
#PBS -M alolex@vcu.edu
#PBS -N STARalign
#PBS -j oe
#PBS -o STARalign-official-081617.log

cd /home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow

INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/input01_toStarAlign.list
OUTDIR=/home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/01_star-align

mkdir -p $OUTDIR

echo "Starting STAR alignment on "`date`

cat $INPUT | while read read1 read2
do
	echo "Starting alignment of $read1 on "`date` 
	prefix=$(basename $read1 .fastq.gz | sed "s/R1_//g") 
	/home/alolex/bin/STAR --runThreadN 48 --genomeDir /home/alolex/refGenomes/STAR_genomes/merged_GRCh38_GRCm38_XMLV --readFilesIn $read1 $read2 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMorder Paired --outReadsUnmapped Fastx --outFileNamePrefix $OUTDIR/$prefix. --quantMode TranscriptomeSAM --outFilterMultimapNmax 1
	echo $prefix
	echo "Finished alignment of $read1 on "`date`

done

echo "Completed STAR alignment on "`date`
