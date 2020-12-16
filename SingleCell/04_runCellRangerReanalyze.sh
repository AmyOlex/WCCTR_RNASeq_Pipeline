#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run the single cell sequencing analysis tool CellRanger Reanalyze from 10X for a list of samples.

REQUIRED ARGUMENTS:
   	-f  Tab delimited input file with list of samples to process and meta data.
	-o	Output directory to save results in.

-f INPUT File:
A tab delimited file with the samples to process along with metadata information and the location of the input raw fastq files.  
Column information must be the following in this order: SampleID, MatrixFile, Paramaters, Barcodes to Keep

-o Directory where the results should be saved. Subfolders will be created using the SampleID given in the input file.


EXAMPLE USAGE:
   >> ~/04_runCellRangerReanalyze.sh -f /path/to/inventoryFile.list -o /path/to/output/dir/

EOF
}

## Parsing input arguments

FILE=
ODIR=
while getopts “hf:o:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
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
if [[ -z $ODIR ]]; then usage; exit 1; fi

cd $ODIR

INPUT=$FILE

cat $INPUT | while read sampleid matrixfile paramfile barcodes
do 
	
	echo cellranger aggr --id=$sampleid --csv=$csvfile --normalize=mapped
	cellranger aggr --id=$sampleid --csv=$csvfile --normalize=mapped

	echo cellranger reanalyze --id=$sampleid --matrix=$matrixfile --params=$paramfile --barcodes=$barcodes
	cellranger reanalyze --id=$sampleid --matrix=$matrixfile --params=$paramfile --barcodes=$barcodes

done
