#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu


## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run Fastq_screen on a list of FastQ files.

REQUIRED ARGUMENTS:
   	-f      Input file with list of FastQ files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 01_cutadapt will be created here.

-f INPUT File:
Input file with list of FastQ files to process.

-d RESULTS Directory:
Directory where results should be saved to.

-c CONFIG File:
Database configuration file.

EXAMPLE USAGE:
   >> ~/00_fastq_screen.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/ -c /path/to/config/file.config

EOF
}

## Parsing input arguments

FILE=
DIR=
CONFIG=
while getopts “hf:d:c:” OPTION
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
	 c)
	     CONFIG=$OPTARG
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
if [[ -z $CONFIG ]]; then usage; exit 1; fi

cd $DIR

INFILE=$FILE

mkdir -p 00_fastqscreen

cat $INFILE | while read file
do
	echo fastq_screen --conf $CONFIG --outdir $DIR/00_fastqscreen/ $file
	fastq_screen --conf $CONFIG --outdir $DIR/00_fastqscreen/ $file
done
