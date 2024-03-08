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
Input File Format | 
Example Usage | ~/02_runCellRanger.sh -f /path/to/sampleFile.list -o /path/to/output/dir/ -n 35 -F
Output Files used in downstream steps? | Yes: PDX Samples have a GEM file used to identify human cells in the step 3, removing mouse cells.  Pure Human/Mouse or Organoid Samples utilize the filtered_bc_matrix data files as input into step 5, finding dead cells.
Output File formats and type | CellRanger “outs” directory with GEM file and "filtered_feature_bc_matrix" directory of files.

## Step 03: Removal of Mouse Cells and Conversion Back to FastQ

Step | Details
---  | ---
Tool Name and Version | subset-bam 1.1.0
Tool Name and Version | bamtofastq v1.3.2
Pipeline Script | 03_create_human_bam.sh
Input File Format | 
Example Usage | ~/03_create_human_bam.sh -f /path/to/sampleFile.list  -d /path/to/output/dir/
Output Files used in downstream steps? | Yes, FastQ files used in step 4, realignment to human genome.
Output File formats and type | BAM and FastQ files that only contain human reads.

## Step 04: Realign to Human only Genome

Step | Details
---  | ---
Tool Name and Version | CellRanger cellranger-6.0.1
Pipeline Script | 02_runCellRanger.sh
Input File Format | See Step 02
Example Usage | See Step 02
Output Files used in downstream steps? | Yes, Human only data: filtered_feature_bc_matrix data files are directly input into step 5, finding dead cells.
Output File formats and type | See Step 02


## Step 05: Finding Dead Cells

Step | Details
---  | ---
Tool Name and Version | R 4.1.3 and hdf5-1.12.2
Pipeline Script | 05_DeadCellAnalysis.R
Input File Format | 
Example Usage | ---NEED TO DOCUMENT---
Output Files used in downstream steps? | Yes, the CSV file “cells2keeps” listing the good quality barcodes is used in step 6, sample merging.
Output File formats and type | Tab delim .txt report, .png images of violin plots before and after filtering, and .csv “cells2keep” files listing barcodes to keep.




































