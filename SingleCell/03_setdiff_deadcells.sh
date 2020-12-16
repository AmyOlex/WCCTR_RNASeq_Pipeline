#!/bin/bash

# Copyright (c) 2019
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to take a list of dead cells and subtract it from a GEM file created during the aligment of the merged human-mouse PDX data.

REQUIRED ARGUMENTS:
	-f	Input metadata file with the following columns for each sample: hg19mm10 sample ID, Path to hg19mm10 outs folder, path to and name of human only dead cell barcode list.

-f Sample list:
Input metadata file with the following columns for each sample: hg19mm10 sample ID, Path to hg19mm10 outs folder, path to and name of human only dead cell barcode list.

EXAMPLE USAGE:
   >> ~/03_create_human_bam.sh -f samples.list

EOF
}

## Parsing input arguments

FILE=
while getopts “hf:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi

cat $FILE | while read sampleid outspath deadpath
do
	echo Processing $sampleid
	echo DEADPATH: $deadpath
	echo GEM: $outspath/analysis/gem_classification.csv
	echo NEWGEM:  $outspath/analysis_deadcells/$sampleid\_gemcellstokeep.csv

	echo COMMAND: grep -F -v -f $deadpath $outspath/analysis/gem\_classification.csv | cut -f 1 -d , > $outspath/analysis_deadcells/$sampleid\_gemcellstokeep.csv
	grep -F -v -f $deadpath $outspath/analysis/gem\_classification.csv | cut -f 1 -d , > $outspath/analysis_deadcells/$sampleid\_gemcellstokeep.csv

done
