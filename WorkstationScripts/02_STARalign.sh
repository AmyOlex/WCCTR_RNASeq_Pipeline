#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run STAR alignment on a list of FastQ files.

REQUIRED ARGUMENTS:
   	-f      Input file with list of FastQ files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 02_STAR will be created here.

-f INPUT File:
Input file with list of FastQ files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/02_STARalign.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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
OUTDIR=$DIR/01_star-align

mkdir -p $OUTDIR

echo "Starting STAR alignment on "`date`

cat $INPUT | while read read1 read2
do
	echo "Starting alignment of $read1 on "`date` 
	prefix=$(basename $read1 .fastq.gz | sed "s/R1_//g") 
	/home/alolex/bin/STAR --runThreadN 48 --genomeDir /home/alolex/refGenomes/STAR_genomes/merged_GRCh38_GRCm38_XMLV --readFilesIn $read1 $read2 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMorder Paired --outReadsUnmapped Fastx --outFileNamePrefix $OUTDIR/$prefix. --quantMode TranscriptomeSAM --outFilterMultimapNmax 1
	echo $prefix
	echo "Finished alignment of $read1 on "`date`

done

echo "Completed STAR alignment on "`date`
