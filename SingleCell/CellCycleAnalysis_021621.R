# Amy Olex
# 2/16/21
# Analysis of Single Cell data to determine Cell Cycle Stage for each cell.
# This Script is to identify the cell cycle stage for each cell and output a csv file for import into Loupe Cell Browser.
#
# This script reads in a config file that has the following information:
# Sample ID
# location of filtered_features_bc_matrix folder for the sample.
# location and file name of the 10X aggr CSV file used to merge the datasets
#
# The output files will be saved to the 10X sudfolder "analysis" with the sampleID as part of the file name.

library("Seurat")
library("readr")
library("png")
library("dplyr")
library("gridExtra")
library("tibble")
library("ggplot2")
library("optparse")

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-c", "--configfile"), type="character", 
              help="Required. Input config file with sample information.", metavar="character"),
  #make_option(c("-m", "--mito"), type="character", 
  #            help="type of mitochondria gene list to use. Options are coding, noncoding, all (default = coding)", 
  #            default = "coding", metavar="character"),
  #make_option(c("-s", "--species"), type="character", 
  #            help="species to use for mitochandria genes. Options are human, mouse, merged (default = human).",
  #            default = "human", metavar="character"),
  #make_option(c("-f", "--mitoFile"), type="character", 
  #            help="The direct path to a mitochondria gene list file.", 
  #            default = NULL, metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character")
  #  make_option(c("--usegem"), type="logical", 
  #              help="flags script to utilize the gem files associated with the inputs to create a second cells2keep file that includes all mouse cells from combined PDX samples.",
  #              default = FALSE, action = "store_true", metavar="logical")
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

## configure the mitochandria file path
#if (is.null(opt$mitoFile)){
#  # If this was not provided, then use the mito and species information.
#  mitoflag = paste0(opt$mito, opt$species)
#  
#  if(mitoflag %in% c("codingmerged", "noncodingmerged")){
#    stop("Mito files for coding and noncoding merged genome files not avaliable.  Use --mitoFile option to specify a file.", call.=FALSE)
#  }
#  mitoFile <- switch(mitoflag, 
#                     "codinghuman" = "/home/scRNASeq/harrell_data/Harrell_SingleCellSequencing/referenceFiles/MitoCodingGenes13_human.txt", 
#                     "noncodinghuman" = "/home/scRNASeq/harrell_data/Harrell_SingleCellSequencing/referenceFiles/MitoNonCodingGenes24_human.txt", 
#                     "allhuman" = "/home/scRNASeq/harrell_data/Harrell_SingleCellSequencing/referenceFiles/MitoAllGenes37_human.txt", 
#                     "codingmouse" = "/home/scRNASeq/harrell_data/Harrell_SingleCellSequencing/referenceFiles/MitoCodingGenes13_mouse.txt", 
#                     "noncodingmouse" = "/home/scRNASeq/harrell_data/Harrell_SingleCellSequencing/referenceFiles/MitoNonCodingGenes24_mouse.txt", 
#                     "allmouse" = "/home/scRNASeq/harrell_data/Harrell_SingleCellSequencing/referenceFiles/MitoAllGenes37_mouse.txt", 
#                     "allmerged" = "/home/scRNASeq/harrell_data/Harrell_SingleCellSequencing/referenceFiles/MitoMasterList_37_hg19mm10.txt")
#  
#  
#} else {
#  mitoFile <- opt$mitoFile
#}

runID <- opt$runid
inFile <- opt$configfile
reportDir <- opt$outdir
reportName <- paste0(reportDir, runID, "_CellCycleStage.csv")

#### Ok, print out all the current running options as a summary.

print("Summary of input options:\n")
print(paste("Run Name:", runID))
print(paste("ConfigFile:", inFile))
print(paste("Report Output:", reportName))
#print(paste("Mitochandria Gene List:", mitoFile))
#print(paste("Using GEM file?", usegem))

# Read in the provided config file and loop for each row.
toProcess = read.table(inFile, header=FALSE, sep="\t")

print(paste(dim(toProcess)[1], " rows were found."))

## Import mito genes
#if(file.exists(mitoFile)){
#  mitogene_ids <- read.delim(mitoFile, header = FALSE, stringsAsFactors = FALSE)[[1]]
#} else {
#  stop(paste0("Mitochondria file provided does not exists: ", mitoFile), call.=FALSE)
#}
#g <- file(reportName, 'w')

#writeLines(c("RunID\tSampleID\tKeptCells\tDeadCells\t%Removed\tMitoCutoff\tlog(nFeatureRange)\tlog(nCountRange)"), g)

for(i in 1:dim(toProcess)[1]){
  print(paste("Processing row", i, "from sample", toProcess[i,1]))
  
  sampleID = as.character(toProcess[i,1])
  datadir = as.character(toProcess[i,2])
  aggrfile = as.character(toProcess[i,3])
  x10dir = paste0(datadir,"filtered_feature_bc_matrix")
  savedir = paste0(datadir, "analysis/")
  # gemfile = paste0(datadir, "analysis/gem_classification.csv")
  
  print(paste("Saving file in: ", savedir))
  #system(paste("mkdir", savedir))
  
  # Open the aggr file
  aggr = read.table(aggrfile, header=TRUE, sep="\t")
  
  # Load the data set and create the Seurat object
  scPDX.data <- Read10X(data.dir = paste0(datadir,"filtered_feature_bc_matrix"))
  # Initialize the Seurat object
  scPDX <- CreateSeuratObject(counts = scPDX.data, project = sampleID, min.cells = 0, min.features = 0)
  
  ########
  ## Rename the cell types based on the CellRanger aggr file
  ########
  
  celltype <- substring(colnames(pbmc@assays$RNA@data), first=18)
  
  print(paste("Testing if aggr library is the same length as unique cell types..."))
  stopifnot(length(aggr$library_id) == length(unique(celltype)))
  print(paste("passed!"))
  
  for(i in 1:length(unique(celltype))){
    celltype[celltype==as.character(i)] <- aggr$library_id[i]
  }
  scPDX$celltype <- celltype
  
  ##########
  ## Run Normalization and Cell Cycle Analysis
  #########
  
  scPDX_regress <- SCTransform(scPDX, vars.to.regress = "nCount_RNA", verbose = FALSE)
  
  scPDX_regress <- CellCycleScoring(scPDX_regress, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes, set.ident = TRUE)
  
  cc <- rownames_to_column(as.data.frame(scPDX_regress$Phase), "Barcode")
  names(cc) <- c("Barcode","CellCyclePhase")
  write.table(x=cc, file=paste0(savedir, runID, "_", sampleID,"_CellCycleStatus.csv"), quote=FALSE, sep=",", row.names = FALSE)
  
}

#close(g)
print("Completed all samples!")










