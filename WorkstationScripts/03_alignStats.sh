#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to Extract alignment statistics from STAR generated BAM files of human/mouse PDX alignments.

REQUIRED ARGUMENTS:
   	-f      Input file with list of PDX BAM files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  The text file output03_bamStats.txt will be created here.

-f INPUT File:
Input file with list of FastQ files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/03_alignStats.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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
OUTFILE=$DIR/output03_bamStats.txt
UMAP1=".Unmapped.out.mate1"
UMAP2=".Unmapped.out.mate2"


touch $OUTFILE

echo "Starting file counting on "`date`
echo -e "Sample,Total Counts,Human Aligned,Mouse Aligned,Unmapped\n" >> $OUTFILE

cat $INPUT | while read file
do
	echo "Starting processing of $file on "`date` 
	dir=`dirname $file`
	base=`basename $file .Aligned.out.bam`
	umapfile1=$dir/$base$UMAP1
	umapfile2=$dir/$base$UMAP2

	
	subtotal=`samtools view $file | wc -l`
	human=`samtools view $file | grep HUMAN | wc -l`
	mouse=`samtools view $file | grep MOUSE | wc -l`
	u1=`cat $umapfile1 | grep @ | wc -l`
	u2=`cat $umapfile2 | grep @ | wc -l`
	total=$(($subtotal+$u1+$u2))
	utotal=$(($u1+$u2))

	echo -e "$base,$total,$human,$mouse,$utotal\n" >> $OUTFILE
	echo "Finished processing of $file on "`date`

done

echo "Completed bam processing on "`date`
