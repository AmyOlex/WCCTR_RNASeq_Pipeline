#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run the single cell sequencing analysis tool CellRanger from 10X for a list of samples.

REQUIRED ARGUMENTS:
   	-f      Tab delimited input file with list of samples to process and meta data.
	-i	Directory where input fastq files are located.
	-o	Output directory to save results in.

-f INPUT File:
A tab delimited file with the samples to process along with metadata information.  Column information must be the following in this order: SampleName, SampleID, Reference, ExpectedCells.
The Reference should be listed as one of the following: hg19, mm10, or hg19mm10

-i Directory where the input fastq files are located.  They should all be int he top level directory and not nested, and must be named using the Illumina default schema.

-o Directory where the results should be saved. Subfolders will be created using the SampleID given in the input file.

EXAMPLE USAGE:
   >> ~/02_runCellRanger.sh -f /path/to/inventoryFile.list -i /path/to/fastq/dir/ -o /path/to/output/dir/

EOF
}

## Parsing input arguments

FILE=
IDIR=
ODIR= 
while getopts “hf:i:o:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
             ;;
         i)
             IDIR=$OPTARG
             ;;
	 o)
	     ODIR=$OPTARG
	     ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi
if [[ -z $IDIR ]]; then usage; exit 1; fi
if [[ -z $ODIR ]]; then usage; exit 1; fi

cd $ODIR

INPUT=$FILE

cat $INPUT | while read sample sampleid ref cellc
do 
	if [[ "$ref" == "mm10" ]]
	then
		reference="/data/refGenomes/CellRanger/refdata-cellranger-mm10-3.0.0"
	elif [[ "$ref" == "hg19mm10" ]]
	then
		reference="/data/refGenomes/CellRanger/refdata-cellranger-hg19-and-mm10-3.0.0"
	else
		usage; exit 1
	fi

	echo Processing $sample with reference $reference;
	cellranger count --id=$sampleid --transcriptome=$reference --fastqs=$IDIR --sample=$sample --expect-cells=$cellc
done
