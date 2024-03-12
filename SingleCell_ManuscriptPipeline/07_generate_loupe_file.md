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
