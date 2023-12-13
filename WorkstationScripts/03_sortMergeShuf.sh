#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run the sort-merge-shuffle process on alignment BAM files that need to be combined.

REQUIRED ARGUMENTS:
   	-f      Input file with list of BAM files to process.

OPTIONAL ARGUMENTS:
	-d	Directory of where results should be saved.  Generally this will be the alignment directory.

-f INPUT File:
Input file with list of BAM files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/03_sortMergeShuf.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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

TMPDIR=~/sort_merge_shuf_tmp_dir

mkdir -p $TMPDIR
mkdir -p $DIR

cd $DIR

echo "Starting sort-merge-shuffle on "`date`
echo "using SAMTools version: " `samtools --version`

cat $FILE | while read file1 file2
do
	echo "Starting sort of $file1 and $file2 on "`date` 
	prefix1=$(basename $file1 .Aligned.toTranscriptome.out.bam) 
	prefix2=$(basename $file2 .Aligned.toTranscriptome.out.bam)

	samtools sort -@ 12 -o $TMPDIR/$prefix1.sort.bam $file1
	samtools sort -@ 12 -o $TMPDIR/$prefix2.sort.bam $file2
	
	samtools merge -r $TMPDIR/$prefix1-$prefix2.merged.bam $TMPDIR/$prefix1.sort.bam $TMPDIR/$prefix2.sort.bam 

	samtools collate $TMPDIR/$prefix1-$prefix2.merged.bam $DIR/$prefix1-$prefix2.merged.shuf.Aligned.toTranscriptome.out

	echo "Finished merging $prefix1 and $prefix2 on "`date`

done

rm -r $TMPDIR

echo "Completed sort-merge-shuffle on "`date`
