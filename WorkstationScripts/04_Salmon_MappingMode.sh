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
OUTDIR=$DIR/04_Salmon_MappingMode

mkdir -p $OUTDIR

echo "Starting Salmon on "`date`

cat $INPUT | while read read1 read2
do
	prefix=$(basename $read1 .R1.nonribo.fq.gc)
	echo "Starting read counting of $prefix on "`date` 
	/usr/local/bin/salmon/bin/salmon quant -i /data/refGenomes/plant/TAIR10.1/TAIR10.1_SalmonIndex -l IU -1 $read1 -2 $read2 --validateMappings -o $OUTDIR/$prefix
	echo $prefix
	echo "Finished read counting of $prefix on "`date`
	mv $OUTDIR/$prefix/quant.sf $OUTDIR/$prefix.quant.sf
done

echo "Completed Salmon on "`date`
