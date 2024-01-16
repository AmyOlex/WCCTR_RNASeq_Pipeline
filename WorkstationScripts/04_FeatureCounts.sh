#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run FeatureCounts on a list of BAM files aligned to the transcriptome.

REQUIRED ARGUMENTS:
   	-f      Input file with list of BAM files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 04_FeatureCounts will be created here.

-f INPUT File:
Input file with list of BAM files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/04_FeatureCounts.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

EOF
}

## Parsing input arguments

FILE=
DIR= 
while getopts “hf:d:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
             ;;
         d)
             DIR=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi
if [[ -z $DIR ]]; then usage; exit 1; fi


INPUT=$FILE
OUTDIR=$DIR/04_FeatureCounts

echo New Directory is $OUTDIR
mkdir -p $OUTDIR

cd $OUTDIR

mamba activate subread

echo "Starting FeatureCounts on "`date`

cat $INPUT | while read file
do
	echo "Starting read counting of $file on "`date` 
	prefix=$(basename $file .Aligned.sortedByCoord.out.bam)
	featureCounts -t gene -g gene_id -p -s 0 --countReadPairs -T 48 -a /vcu_gpfs2/home/harrell_lab/refGenomes/merged_GDC_GRCh38_GRCm38_XMLV/merged_gencode_humanV22_mouseVM12.annotation.gtf -o $OUTDIR/$file_featureCounts.txt $file
	
	echo $prefix
	echo "Finished read counting of $file on "`date`
done

echo "Completed Salmon on "`date`
