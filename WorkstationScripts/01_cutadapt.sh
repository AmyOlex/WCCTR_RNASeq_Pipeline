#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run CutAdapt on a list of FastQ files.

REQUIRED ARGUMENTS:
   	-f      Input file with list of FastQ files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 01_cutadapt will be created here.

-f INPUT File:
Input file with list of FastQ files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/01_cutadapt.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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

INFILE=$FILE

mkdir -p 01_cutadapt

cat $INFILE | while read s1 s2
do
	base1=`basename $s1 .fastq.gz`
	base2=`basename $s2 .fastq.gz`
	
	## Illumina TruSeq Universal Adaptors
	#cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -u 5 -U 5 -o 01_cutadapt/$base1.trimmed.fastq.gz -p 01_cutadapt/$base2.trimmed.fastq.gz $s1 $s2 --minimum-length 30

	## Illumina Universal Adaptors
	cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o 01_cutadapt/$base1.trimmed.fastq.gz -p 01_cutadapt/$base2.trimmed.fastq.gz -j 30 $s1 $s2 --minimum-length 30	




done


#Illumina Universal Adapter: AGATCGGAAGAG 
#Illumina Small RNA 3' Adapter: TGGAATTCTCGG 
#Illumina Small RNA 5' Adapter: GATCGTCGGACT 
#Nextera Transposase Sequence: CTGTCTCTTATA
#SOLID Small RNA Adapter: CGCCTTGGCCGT
