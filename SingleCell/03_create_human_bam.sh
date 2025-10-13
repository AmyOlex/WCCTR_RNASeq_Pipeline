#!/bin/bash

# Copyright (c) 2019
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used on the Cell Ranger count output and must be run from the outs directory.  It extracts all the human barcodes from teh GEM analysis, filters the BAM file, sorts the BAM file and then exports to 2 Fastq files.

REQUIRED ARGUMENTS:
	-f	Input metadata file with the following columns for each sample: SampleName, PathToGEM, PathToBAM
OPTIONAL ARGUMENTS:
	-d		Directory of where the resulting fastq files should be saved.  Default is the current directory.
	-G      The alias to use for identifying Graft (human) barcodes.  Default is GRCh38.
	-H      The alias to use for identifying Host (mouse) barcodes.  Default is mm10.
	-M      The alias to use for identifying Multiplet barcodes.  Default is Multiplet.

-f Sample list:
Tab delimited list of samples to process. Each row includes a unique sample name followed by the path to the GEM file for that sample and the the path to the BAM file for that sample.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/03_create_human_bam.sh -f samples.list  -d /path/to/results/dir/

EOF
}

## Parsing input arguments

FILE=
DIR=
GRAFT=
HOST=
MULTI=
while getopts “hf:d:G:H:M:” OPTION
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
         G)
             GRAFT=$OPTARG
             ;;
         H)
             HOST=$OPTARG
             ;;
         M)
             MULTI=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi
if [[ -z $DIR ]]; then DIR="."; fi
if [[ -z $GRAFT ]]; then GRAFT="GRCh38"; fi
if [[ -z $HOST ]]; then HOST="mm10"; fi
if [[ -z $MULTI ]]; then MULTI="Multiplet"; fi

echo "Using Graft Alias: $GRAFT"
echo "Using Host Alias: $HOST"
echo "Using Multiplet Alias: $MULTI"

cat $FILE | while read sampleid gempath bampath
do
	echo Filtering Barcodes for $sampleid using GEM at $gempath

	tail -n +2 $gempath/gem_classification.csv | grep $GRAFT | cut -d ',' -f1 > $gempath/graft_barcodes.txt 
	tail -n +2 $gempath/gem_classification.csv | grep $HOST | cut -d ',' -f1 > $gempath/host_barcodes.txt
	tail -n +2 $gempath/gem_classification.csv | grep $MULTI | cut -d ',' -f1 > $gempath/multiplet_barcodes.txt

	echo Subsetting BAM

	/lustre/home/harrell_lab/src/subset-bam_linux --bam $bampath/possorted_genome_bam.bam --cell-barcodes $gempath/graft_barcodes.txt --out-bam $bampath/graft_genome_bam.bam --cores 48

	echo Converting BAM to Fastq

	/lustre/home/harrell_lab/src/bamtofastq-1.3.2 $bampath/graft_genome_bam.bam $DIR/$sampleid

	echo Renaming and Moving Files

	cd $DIR/$sampleid
	mv */*.gz .
	#newdir=`ls`
	#cd $newdir
	rename bamtofastq $sampleid bamtofastq*

	##mv * $DIR/
	echo Completed $sampleid
done
