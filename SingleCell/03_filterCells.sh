#!/bin/bash

# Copyright (c) 2019
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is written to use a list of cell barcodes to keep and extracts them from the sample BAM file, then converts to fastq.

REQUIRED ARGUMENTS:
	-f	Input metadata file with the following columns for each sample: SampleName, PathToCellsToKeepFile, PathToBAM
OPTIONAL ARGUMENTS:
	-d	Directory of where the resulting fastq files should be saved.  Default is the current directory.

-f Sample list:
Tab delimited list of samples to process. Each row includes a unique sample name followed by the path to the cells to keep for this sample followed by the path to the BAM file for that sample.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/03_filterCells.sh -f samples.list  -d /path/to/results/dir/

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
if [[ -z $DIR ]]; then DIR="."; fi

cat $FILE | while read sampleid cellpath bampath
do
	echo Filtering $sampleid using barcodes from file $cellpath

	echo Subsetting BAM

	subset-bam --bam $bampath/possorted_genome_bam.bam --cell-barcodes $cellpath --out-bam $bampath/$sampleid_filtered.bam --cores 20

	echo Converting BAM to Fastq

	bamtofastq $bampath/$sampleid_filtered.bam $DIR/$sampleid

	echo Renaming and Moving Files

	cd $DIR/$sampleid
	newdir=`ls`
	cd $newdir
	rename bamtofastq $sampleid bamtofastq*

	##mv * $DIR/
	echo Completed $sampleid
done
