#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M alolex@vcu.edu
#PBS -N Salmon2
#PBS -j oe
#PBS -o Salmon2.log


INPUT=/home/alolex/data/Harrell/WHIM2_altPipeline_8-23-17/whim2_mouse_bam_for_salmon.list
OUTDIR=/home/alolex/data/Harrell/WHIM2_altPipeline_8-23-17/salmon_output_mouse

mkdir -p $OUTDIR

echo "Starting Salmon on "`date`

cat $INPUT | while read file
do
	echo "Starting read counting of $file on "`date` 
	prefix=$(basename $file .Aligned.toTranscriptome.out.mouse.bam)
	salmon quant -t ~/refGenomes/merged_GDC_GRCh38_GRCm38_XMLV.transcripts.fa -l IU -a $file -o $OUTDIR/$prefix
	echo $prefix
	echo "Finished read counting of $file on "`date`
	mv $OUTDIR/$prefix/quant.sf $OUTDIR/$prefix/$prefix.mouse.quant.sf
done

echo "Completed Salmon on "`date`
