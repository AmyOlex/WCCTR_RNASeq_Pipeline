# Single Cell RNASeq Pipeline used in Altman Manuscript
Current Pipeline used to process and integrate the latest 2023 samples and organoids into the master merge and collaboration Visvader master merge with subtyping analysis.

## Step 01: FastQC and MultiFastQC

Step | Details
---  | ---
Tool Name and Version | FastQC v0.11.9
Pipeline Script | 01_fastqc.sh
Input File Format | Text file with one FastQ file name and absolut path per line.
Example Usage | ~/01_fastqc.sh -f /path/to/fastqFiles.list -d /path/to/output/dir/

Step | Details
---  | ---
Tool Name and Version | MultiQC, version 1.17
Pipeline Script | None, manually run MultiQC in same directory as FastQC output.
Output Files used in downstream steps? | No, FastQC and MultiQC reports are manually reviewed to inform on data quality.


## Step 02: Single Sample Alignment

Step | Details
---  | ---
Tool Name and Version | CellRanger, cellranger-6.0.1
Pipeline Script |  02_runCellRanger.sh
Input File Format | No header, TAB delimited, has the following 5 columns: sample ID used in fastq file name, result directory name, reference genome (grch38, mm10, or grch38mm10), number of cells to target, absolute path to raw fastq data directory.
Example Usage | ~/02_runCellRanger.sh -f /path/to/sampleFile.list -o /path/to/output/dir/ -n 35
Output Files used in downstream steps? | Yes: PDX Samples have a GEM file used to identify human cells in the step 3, removing mouse cells.  Pure Human/Mouse or Organoid Samples utilize the filtered_bc_matrix data files as input into step 5, finding dead cells.
Output File formats and type | CellRanger “outs” directory with GEM file and "filtered_feature_bc_matrix" directory of files.

## Step 03: Removal of Mouse Cells and Conversion Back to FastQ

Step | Details
---  | ---
Tool Name and Version | subset-bam 1.1.0
Tool Name and Version | bamtofastq v1.3.2
Pipeline Script | 03_create_human_bam.sh
Input File Format | TAB delimited, has the following 3 columns: sample ID (same as Step 02), absolute path to samples “outs/analysis” directory, absolute path to samples “outs” directory.
Example Usage | ~/03_create_human_bam.sh -f /path/to/sampleFile.list  -d /path/to/output/dir/
Output Files used in downstream steps? | Yes, FastQ files used in step 4, realignment to human genome.
Output File formats and type | BAM and FastQ files that only contain human reads.

## Step 04: Realign to Human only Genome

Step | Details
---  | ---
Tool Name and Version | CellRanger cellranger-6.0.1
Pipeline Script | 02_runCellRanger.sh
Input File Format | See Step 02; must edit the number of cells to those that are human only instead of the default targeted number.  Also add the "-F" for force cells option to stop CellRanger from doing further filtering of cells.
Example Usage | ~/02_runCellRanger.sh -f /path/to/sampleFile.list -o /path/to/output/dir/ -n 35 -F
Output Files used in downstream steps? | Yes, Human only data: filtered_feature_bc_matrix data files are directly input into step 5, finding dead cells.
Output File formats and type | See Step 02


## Step 05: Finding Dead Cells

Step | Details
---  | ---
Tool Name and Version | R 4.1.3 and hdf5-1.12.2
Pipeline Script | 05_DeadCellAnalysis.R
Input File Format | TAB delimited, has the following 2 columns: sample ID, absolute path to human-aligned “outs” directory.
Example Usage (use --help option for all parameters) | Rscript 05_DeadCellAnalysis.R -r UniqRunID -c sample_tabdelim_file.list -f /path/to/mitogenes/MitoCodingGenes13_human.txt -s human -o /path/to/results/directory
Output Files used in downstream steps? | Yes, the CSV file “cells2keeps” listing the good quality barcodes is used in step 6, sample merging.
Output File formats and type | Tab delim .txt report, .png images of violin plots before and after filtering, and .csv “cells2keep” files listing barcodes to keep.

## Step 06: Sample Merging

Step | Details
---  | ---
Tool Name and Version | R 4.1.3 and hdf5-1.12.2
Pipeline script | 06_SeuratMerge.R
Input File Format | Comma delimited, has header row (SampleName,DataType,SamplePath,Source,Condition,Sex,TumorType,Cells2Keep), has the following 8 columns (but more can be added and will be used as cell annotations): Sample ID, data type (10X or seurat), absolute path to sample’s “filtered_feature_bc_matrix” directory, sample source (e.g. PDX, MGT, etc), condition or treatment, sample sex (M or F), tumor type (e.g. HCI011, UCD52, etc), absolute path to CSV file listing barcodes of cells to keep in merge.
Example Usage (use --help option for all parameters) | Rscript 06_SeuratMerge.R -r UniqRunID -c /path/to/input.csv -i LogNormalize -t simple --parallel -n 32 --filter
Output Files used in downstream steps? | Yes, .h5 file is used by CellRanger Reanalyze to generate a Loupe file
Output File formats and type | .h5 file, RData file contining the merged Seurat object, multiple .CSV annotation files (one for each column in the input CSV file)

## Step 07: Create Loupe File

Note: During the development of this manuscript the RLoupe package that converts Seurat objects to Loupe-compatible files had not been released and this was our work-around. 

Step | Details
---  | ---
Tool Name and Version | CellRanger cellranger-6.0.1
Pipeline script | N/A, calls the reanalyze function directly.
Input File Format | the .h5 file generated from the merging script in step 06.
Example Usage | cellranger reanalyze --id=UniquRunID --description=RunDescription --matrix=/path/to/merged_file.h5 --params=/path/to/reanalyze_params.txt
Output Files used in downstream steps? | Yes, Loupe file used for figure generation and analysis.
Output File formats and type | CellRanger “outs” directory with a Loupe File.








































