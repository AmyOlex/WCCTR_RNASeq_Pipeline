#!/bin/bash

# Copyright (c) 2019
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run BBduk.sh on a list of FastQ files to remove human rRNA contaminated reads.

REQUIRED ARGUMENTS:
   	-f      Input file with list of FastQ files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 01_rrna_filtered will be created here.

-f INPUT File:
Input file with list of FastQ files to process.

-d RESULTS Directory:
Directory where results should be saved to.

-r REFERENCE contaminant file.
Path to and name of file containing the rRNA sequences.

EXAMPLE USAGE:
   >> ~/01_rRNAFilter.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/ -r /path/to/rRNA.fasta

EOF
}

## Parsing input arguments

FILE=
DIR = 
REF = 
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

mkdir -p 01_rrna_filtered

OUTDIR=$DIR/01_rrna_filtered

cat $INPUT | while read file
do 
	prefix=$(basename $file .fastq.gz)
	echo "Processing $file ...";
	bbduk.sh in=$file outm=$OUTDIR/$prefix.ribo.fq outu=$OUTDIR/$prefix.nonribo.fq ref=$REF stats=$OUTDIR/$prefix.rRNAStats.txt
	gzip $OUTDIR/$prefix.ribo.fq
	gzip $OUTDIR/$prefix.nonribo.fq
done
