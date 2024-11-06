## Amy Olex
## 09/20/2024
## Converts a mouse Seurat object to human and then saves as a h5ad file for import into Python.
## This is data re-formatting for the LIANA Cell-Cell Communication analysis.

#### Required Software Installation

#BiocManager::install("rhdf5")
#BiocManager::install("ensembldb")
#devtools::install_github("AmyOlex/seuratTools", ref="debug_cross_species")
#devtools::install_github("AmyOlex/SeuratHelper", ref="debug_write_h5ad")
#install.packages("XML")
#BiocManager::install("wiggleplotr")
#BiocManager::install("velociraptor")

## Load packages
#library(velociraptor)
#library(wiggleplotr)
#library(GenomicFeatures)
#library(rtracklayer)
#library(biomaRt)
#library(restfulr)
#library(ensembldb)

library(rhdf5)
library(seuratTools)
library(SeuratHelper)
library(Seurat)
library(optparse)

## Tested successfully on the following dataset:
#runID<- "SimpleMerge_Mouse_UntreatedOnly_noForce_Log2_GRCh38_240830_Seurat_simpleMerge_LogNorm"
#outDir <- "/lustre/home/alolex/SwatiDeb_R01_Cell2Cell-Communication/data/"
#inFile <- "/lustre/home/alolex/SwatiDeb_R01_Cell2Cell-Communication/data/SimpleMerge_Mouse_noForce_Log2_GRCh38_UntreatedOnly_240830_Seurat_simpleMerge_LogNormalize_Annotated.rds"
#slot = "data"

options(future.globals.maxSize = 100000 * 1024^2)

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-i", "--input"), type="character", 
              help="Required. Input RDS file with the saved Seurat object.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-s", "--slot"), type="character", 
              help="The name of the data slot to convert to h5ad format. Default is counts.",
              default = "counts", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("An input RDS file with a Seurat object must be provided.", call.=FALSE)
}

runID <- opt$runid
inFile <- opt$input
outDir <- opt$outdir
slot <- opt$slot


## load data
seurat.obj <- readRDS(inFile)
seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize")

## Convert Mouse to Human
seurat.obj_human <- seuratTools::convert_mouse_seu_to_human(seurat.obj)

## I may implement this code in the future, but right now it is not recommended 
##. because it does not guarantee all cell types will be sampled, which causes issues with Cell-Cell Communication packages.
## Subset the Seurat Object for testing
#samplesize <- floor(0.1*length(names(seurat.obj_human$orig.ident)))
#set.seed(42)
#sampledcells <- sample(x = names(seurat.obj_human$orig.ident), size = samplesize, replace = F)
#seurat.sub <- SeuratObject:::subset.Seurat(seurat.obj_human, cells = sampledcells)

## Write the Seurat object to h5ad files
SeuratHelper::write_h5ad(object = seurat.obj_human, slot=slot, file = paste0(outDir, runID, "_Mouse2HumanGeneID.h5ad"))
SeuratHelper::write_h5ad(object = seurat.obj, slot=slot, file = paste0(outDir, runID, "_MouseGeneID.h5ad"))




