## Amy Olex
## 01.16.24
## featureCounts pipeline for merging read counts from BulkRNASeq data with multiple samples and calculating the TPM values.  
##
## This script imports the human, mouse, or pdx data for all samples from the featureCount files. 
## It extracts the read counts and recalculates the TPM values for human and mouse separately.
## Gene names and the preferred sample names are appended to the TPM and raw count data frames.
## Three data files are then written out for each species: gene metadata, count data, TPM values.
## This script does nothing else except merge and format the raw expression values and calculate the TPM values for mouse and human individually.

#library(tximport)
library(readr)
library(AnnotationHub)
#library(ensembldb)
#library(RNASeqBits)
library(NMF)
#library(limma)
library(optparse)

calc.tpm.fromFeatureCounts <- function(metadata, data){
  
  RPK <- matrix(0, nrow=dim(data)[1], ncol=dim(data)[2])
  
  for(row in 1:dim(data)[1]){
    for(col in 1:dim(data)[2]){
      RPK[row,col] <- data[row,col]/metadata$Length[row]
    }
  }
  
  ##Calculate the sums of each column and divide by 1000000
  scale_factor <- colSums(RPK)/1000000
  
  ##Now divide all values in each column by the scaling factor
  TPM <- t(t(RPK)/scale_factor)
  colnames(TPM) <- colnames(data)
  row.names(TPM) <- rownames(data)
  return(as.data.frame(TPM))
}

localtest = FALSE
###########################################
#### Local Testing Block
if(localtest){
  setwd("/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/featureCount_testing")
  runID <- "TEST"
  inFile <- "./merge_counts_testing_config.txt"
  outDir <- "./"
  species <- "pdx"
  savedir <- paste0(outDir,runID)
  IdCol <- "Geneid"
  options(future.globals.maxSize = 3000 * 1024^2)
}
###########################################


option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this dataset.", metavar="character"),
  make_option(c("-c", "--configfile"), type="character", 
              help="Required. Input config file with sample information. Should include the columns: Name, File", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="Output directory to save results (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-s", "--species"), type="character", 
              help="The data species. If PDX then both human and mouse will be utilized (HUMAN, MOUSE, PDX). Default is human.",
              default = "HUMAN", metavar="character"),
  make_option(c("-i", "--idcol"), type="character", 
              help="The name of the ID column that contains Ensemble IDs. Default = Geneid",
              default = "Geneid", metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$configfile)){
  print_help(opt_parser)
  stop("A configuration file with sample information must be provided.", call.=FALSE)
}


runID <- opt$runid
inFile <- opt$configfile
outDir <- opt$outdir
savedir <- paste0(outDir,runID)
species <- opt$species
IdCol <- opt$idcol

#### Print out all the current running options as a summary.

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("ConfigFile: ", inFile))
print(paste("OutDir: ", outDir))
print(paste("Species: ", species))
print(paste("SaveDir: ", savedir))
print(paste("IdCol: ", IdCol))

############## PROCESS CONFIG FILE###########################
# Read in the provided config file and loop for each row.
toProcess = read.table(inFile, header=TRUE, sep=",", stringsAsFactors = FALSE)

## check to ensure they exist
stopifnot(all(file.exists(toProcess$File)))

## initiate the merged count variable
merged_counts <- ""

## Import Annotation DB to get gene Symbol identifiers.
ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")

if(species %in% c("PDX", "HUMAN")){ orgdb_h <- query(orgs, "Homo sapiens")[[1]] }
if(species %in% c("PDX", "MOUSE")){ orgdb_m <- query(orgs, "Mus musculus")[[1]] }

## Process each file and merge it's read counts into the same dataframe
for (i in seq_along(toProcess$File)) {
  message(i, " ", appendLF = FALSE)
  raw <- read.table(toProcess$File[i], header=TRUE, sep="\t", stringsAsFactors = FALSE)
  names(raw)[7] <- make.names(toProcess$Name[i])
  
  if (i == 1) {
    ## get all the annotations I will need imported.
    merged_counts <- raw[,1:6]
    ## Remove the version identifier from the Ensemble ID
    merged_counts$Ensembl = gsub("\\..*$", "",merged_counts[[IdCol]])
    
    geneSym_h <- NULL
    geneSym_m <- NULL
    
    if(species=="PDX"){
      merged_counts$species <- sapply(merged_counts$Chr, function(x) strsplit(x, '_')[[1]][1])
    }
    else{
      merged_counts$species <- rep(species,nrow(merged_counts))
    }

    if(species %in% c("PDX", "HUMAN")){
      k_h <- merged_counts[merged_counts$species == "HUMAN","Ensembl"]
      geneSym_h <- as.data.frame(mapIds(orgdb_h, keys=k_h, column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
      names(geneSym_h) <- c("geneSym")
    }
    
    if(species %in% c("PDX", "MOUSE")){
      k_m <- merged_counts[merged_counts$species == "MOUSE","Ensembl"]
      geneSym_m <- as.data.frame(mapIds(orgdb_m, keys=k_m, column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
      names(geneSym_m) <- c("geneSym")
    }
    
    geneSym <- rbind(geneSym_h, geneSym_m)
    stopifnot(all(merged_counts$Ensembl == row.names(geneSym)))
    
    merged_counts$geneSym <- geneSym[,1]
    
  }
  stopifnot(all(merged_counts[[IdCol]] == raw[[IdCol]]))
  merged_counts <- merge(x = merged_counts, y = raw[,c(1,7)], by = IdCol, all.x=TRUE, sort=FALSE)
  
}

# Assign Row Names
row.names(merged_counts) <- merged_counts$Ensembl


if(species %in% c("PDX", "HUMAN")){
  # Extract metadata columns and count columns
  merged_meta_h <- merged_counts[merged_counts$species == "HUMAN", 1:9]
  merged_counts_h <- merged_counts[merged_counts$species == "HUMAN",10:length(merged_counts)]
  merged_counts_h_annot <- cbind(merged_meta_h[,"Ensembl",drop=FALSE], merged_counts_h)
  
  # Calculate the TPM values
  TPM_h <- calc.tpm.fromFeatureCounts(merged_meta_h, merged_counts_h)
  TMP_h_annot <- merge(merged_meta_h[,"Ensembl", drop=FALSE], TPM_h, by="row.names")[,-1]
  
  # Save files
  write.table(TMP_h_annot, file=paste0(savedir, "_Human_TPM.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  write.table(merged_counts_h_annot, file=paste0(savedir, "_Human_counts.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  write.table(merged_meta_h, file=paste0(savedir, "_Human_geneMetaData.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  
  
}

if(species %in% c("PDX", "MOUSE")){
  # Extract metadata columns and count columns
  merged_meta_m <- merged_counts[merged_counts$species == "MOUSE", 1:9]
  merged_counts_m <- merged_counts[merged_counts$species == "MOUSE",10:length(merged_counts)]
  merged_counts_m_annot <- cbind(merged_meta_m[,"Ensembl",drop=FALSE], merged_counts_m)
  
  # Calculate the TPM values
  TPM_m <- calc.tpm.fromFeatureCounts(merged_meta_m, merged_counts_m)
  TMP_m_annot <- merge(merged_meta_m[,"Ensembl", drop=FALSE], TPM_m, by="row.names")[,-1]
  
  # Save files
  write.table(TMP_m_annot, file=paste0(savedir, "_Mouse_TPM.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  write.table(merged_counts_m_annot, file=paste0(savedir, "_Mouse_counts.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  write.table(merged_meta_m, file=paste0(savedir, "_Mouse_geneMetaData.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  
}

print("COMPLETED")

