## Amy Olex
## 12.18.24
## Pipeline for merging read counts from BulkRNASeq results obtained from the JAX PDXNet NextFlow output.  
##
## This script imports the human, mouse, or pdx data for all samples from the RSEM "genes.results" files.
## Gene names and the preferred sample names are appended to the TPM and raw count data frames.
## Three data files are then written out for each species: gene metadata, count data, TPM values.
## This script does nothing else except merge and format the raw expression values and calculate the TPM values for mouse and human individually.

#library(readr)
#library(AnnotationHub)
#library(NMF)
library(optparse)
#library(dbplyr)
library(tools)

sessionInfo()

## Initializes the counts matrix to have all the annotation columns from Feature Counts.
init.merged.counts <- function(fname){
  
  ## load in the first file of the list
  raw <- read.table(fname, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  
  ## get all the annotations columns.
  merged_counts <- raw[,1:4]
  
  ## parse out gene symbols
  merged_counts$geneSym <- sapply(raw$gene_id, function(x) strsplit(x, '_')[[1]][2])
  
  return(merged_counts)
}

localtest = FALSE
###########################################
#### Local Testing Block
if(localtest){
  setwd("~/test/")
  runID <- "TEST"
  inDir <- "/lustre/home/harrell_lab/bulkRNASeq/jax_pipeline_results/human_isoform_exp/"
  outDir <- "~/test/"
  species <- "PDX"
  savedir <- paste0(outDir,runID)
  #IdCol <- "gene_id"
  IdCol <- "transcript_id"
  ctype <- "ISOFORM"
  options(future.globals.maxSize = 3000 * 1024^2)
}
###########################################


option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this dataset.", metavar="character"),
  #make_option(c("-c", "--configfile"), type="character", 
  #            help="Required. Input config file with sample information. Should include the columns: Name, File", metavar="character"),
  make_option(c("-d", "--indir"), type="character", 
              help="Required. Input directory with files, or symlinks to files to include in the merged output.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="Output directory to save results (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-s", "--species"), type="character", 
              help="The data species. If PDX then both human and mouse will be utilized (HUMAN, MOUSE, PDX). Default is human.",
              default = "HUMAN", metavar="character"),
  make_option(c("-t", "--counttype"), type="character", 
              help="The type of count data to generate (GENE, ISOFORM). Isoform is only valid for RSEM count data output. Default is GENE",
              default = "GENE", metavar="character"),
  make_option(c("-i", "--idcol"), type="character", 
              help="The name of the ID column that contains unique row names. For Salmon quant data it should be Geneid, for RSEM counts use the default gene_id for gene expression, or use transcript_id for isoform counts. Default = gene_id",
              default = "gene_id", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$indir)){
  print_help(opt_parser)
  stop("An input directory where count files are located is required.", call.=FALSE)
}


runID <- opt$runid
inDir <- opt$indir
outDir <- opt$outdir
savedir <- paste0(outDir,"/",runID)
species <- opt$species
IdCol <- opt$idcol
ctype <- opt$counttype

#### Print out all the current running options as a summary.

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("DataDir: ", inDir))
print(paste("OutDir: ", outDir))
print(paste("Species: ", species))
print(paste("SaveDir: ", savedir))
print(paste("IdCol: ", IdCol))
print(paste("CountType: ", ctype))

############## PROCESS CONFIG FILE###########################
# Read in the provided config file and loop for each row.
toProcess = data.frame(File = file.path(inDir, list.files(inDir))) #read.table(inFile, header=TRUE, sep=",", stringsAsFactors = FALSE)

## check to ensure they exist
stopifnot(all(file.exists(toProcess$File)))
toProcess$Name <- sub("_(human|mouse)\\.(genes|isoforms)\\.results$", "", basename(toProcess$File))



## initiate the merged count matrix with the gene annotations
merged_counts <- init.merged.counts(toProcess$File[1])
merged_tpm <- merged_counts
merged_fpkm <- merged_counts
merged_iso <- merged_counts

print("Processing Samples:")

for(i in seq_along(toProcess$File)){
  message(toProcess$Name[i], " ", appendLF = FALSE)
  print(toProcess$Name[i])
  rawdata <- read.table(toProcess$File[i], header = TRUE)
  row.names(rawdata) <- rawdata[,IdCol]
  
  countdata <- rawdata[,"expected_count", drop=FALSE]
  names(countdata) <- toProcess$Name[i]
  tpmdata <- rawdata[,"TPM", drop=FALSE]
  names(tpmdata) <- toProcess$Name[i]
  fpkmdata <- rawdata[,"FPKM", drop=FALSE]
  names(fpkmdata) <- toProcess$Name[i]
  
  merged_counts <- merge(x = merged_counts, y = countdata, by.x = IdCol, by.y = "row.names", all.x=TRUE, sort=FALSE)
  merged_tpm <- merge(x = merged_tpm, y = tpmdata, by.x = IdCol, by.y = "row.names", all.x=TRUE, sort=FALSE)
  merged_fpkm <- merge(x = merged_fpkm, y = fpkmdata, by.x = IdCol, by.y = "row.names", all.x=TRUE, sort=FALSE)
  
  if(ctype=="ISOFORM"){
    isodata <- rawdata[,"IsoPct", drop=FALSE]
    names(isodata) <- toProcess$Name[i]
    merged_iso <- merge(x = merged_iso, y = isodata, by.x = IdCol, by.y = "row.names", all.x=TRUE, sort=FALSE)
  }
}

# Save files
write.table(merged_counts, file=paste0(savedir, "_counts.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(merged_tpm, file=paste0(savedir, "_TPM.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(merged_fpkm, file=paste0(savedir, "_FPKM.tsv"), quote=FALSE, sep="\t", row.names=FALSE)

if(ctype=="ISOFORM"){
  write.table(merged_iso, file=paste0(savedir, "_IsoPct.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  
}

print(paste("COMPLETED ", runID))
print(paste("Number of samples processed: ", length(toProcess$File)))


