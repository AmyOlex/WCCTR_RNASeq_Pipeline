#!/bin/bash

# Copyright (c) 2021
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run SAMTools-sort on a list of BAM files to sort them by coordinates.

REQUIRED ARGUMENTS:
   	-f      Input file with list of BAM files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 03_coordSorted will be created here.

-f INPUT File:
Input file with list of BAM files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/03_coordSort.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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


cd $DIR

INPUT=$FILE

mkdir -p 03_coordSorted

OUTDIR=$DIR/03_coordSorted

cat $INPUT | while read file
do 
	echo $file;
	samtools sort -o $file.sorted.bam -O bam -@ 30 $file
done
