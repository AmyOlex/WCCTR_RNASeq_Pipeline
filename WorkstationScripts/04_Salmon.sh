#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run Salmon on a list of BAM files aligned to the transcriptome.

REQUIRED ARGUMENTS:
   	-f      Input file with list of BAM files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 04_Salmon will be created here.

-f INPUT File:
Input file with list of BAM files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/04_Salmon.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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
OUTDIR=$DIR

mkdir -p $OUTDIR

cd $OUTDIR

echo "Starting Salmon on "`date`

cat $INPUT | while read file1 file2
do
	echo "Starting read counting of $file1 on "`date` 
	prefix=$(basename $file1 .Aligned.toTranscriptome.out.bam)
	/vcu_gpfs2/home/mccbnfolab/src/salmon-1.5.2_linux_x86_64/bin/salmon quant -t /vcu_gpfs2/home/mccbnfolab/ref_genomes/refgenie/alias/hg38/salmon_sa_index/default/ -l A -1 $file1 -2 $file2 -o $OUTDIR/$prefix
	echo $prefix
	echo "Finished read counting of $file1 on "`date`
done

echo "Completed Salmon on "`date`
