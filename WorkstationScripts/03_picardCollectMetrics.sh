#!/bin/bash

# Copyright (c) 2021
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run SAMTools-stats on a list of BAM files.

REQUIRED ARGUMENTS:
   	-f      Input file with list of BAM files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 03_picardAlignmentMetrics will be created here.

-f INPUT File:
Input file with list of BAM files to process.

-d RESULTS Directory:
Directory where results should be saved to.

-r REFERENCE Fasta:
Path to the reference fasta file used for alignment.

EXAMPLE USAGE:
   >> ~/03_picardCollectMetrics.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/ -r /path/to/reference.fasta

EOF
}

## Parsing input arguments

FILE=
DIR=
REF=
while getopts “hf:d:r:” OPTION
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
         r)
             REF=$OPTARG
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
if [[ -z $REF ]]; then usage; exit 1; fi

cd $DIR

INPUT=$FILE

mkdir -p 03_picardAlignmentMetrics

OUTDIR=$DIR/03_picardAlignmentMetrics

cat $INPUT | while read file
do 
	echo $file;
	java -jar /vcu_gpfs2/home/mccbnfolab/src/picard.jar CollectAlignmentSummaryMetrics I=$file O=$file.stats R=$REF
done
