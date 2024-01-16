## scRNASeq Annotation Porting
## Transfer annotations from a superset to a subset or vice versa.
## Amy Olex 
## 9/14/23

library("optparse")
library("dplyr")
library("reshape2")


############
## Helper functions
###########

options(future.globals.maxSize = 100000 * 1024^2)

start_time <- Sys.time()

#infile = "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Pseudobulk_PAM50_subtyping/100422_Visvader-MiniSubset_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
#WD = "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Pseudobulk_PAM50_subtyping/"
#setwd(WD)

option_list = list(
  make_option(c("-s", "--sourcelib"), type="character", default=NULL, 
              help="Required. The name and path to the Source dataset LibraryID CSV annotation file. 2-column format with barcode and LibraryID.", metavar="character"),
  make_option(c("-t", "--targetlib"), type="character", 
              help="The name and path to the Target dataset LibraryID CSV annotation file. 2-column format with barcode and LibraryID.", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", 
              help="Annotation that needs to be ported. 2-column format with barcode and LibraryID.",
              default = "orig.ident", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-f", "--outfile"), type="character", 
              help="output file name for results report. Default is ported.csv.",
              default = "ported.csv", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$sourcelib)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$targetlib)){
  print_help(opt_parser)
  stop("An input RData file with Seurat single cell object must be provided.", call.=FALSE)
}

## Cancer Only 
#sourcelib <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#outDir <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/"

## to TN Cancer Only - pseudo PAM50
#outFile <- "050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_pseudoPAM50_ported-to-TNCancerOnly.csv"
#targetlib <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/TN_051023_Visvader-plus-PDX_MASTER_CancerTNOnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#annot <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_pseudoPAM50_Annotation.csv"

## to TN Cancer Only - pseudo Claudin Low
#outFile <- "050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_ClaudinLow_ported-to-TNCancerOnly.csv"
#targetlib <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/TN_051023_Visvader-plus-PDX_MASTER_CancerTNOnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#annot <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_pseudoClaudinLow_Annotation.csv"

## to TN Cancer Only - SCSubtype
#outFile <- "050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_SCSubtype_ported-to-TNCancerOnly.csv"
#targetlib <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/TN_051023_Visvader-plus-PDX_MASTER_CancerTNOnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#annot <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/SCSubtyper_testing/062823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_SCSubtypes_Annotation.csv"




## to ER Cancer Only - pseudo PAM50
#outFile <- "050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_pseudoPAM50_ported-to-ERCancerOnly.csv"
#targetlib <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/CancerOnlySubsets/ER_051023_Visvader-plus-PDX_MASTER_CancerEROnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#annot <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_pseudoPAM50_Annotation.csv"

## to ER Cancer Only - pseudo Claudin Low
#outFile <- "050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_ClaudinLow_ported-to-ERCancerOnly.csv"
#targetlib <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/CancerOnlySubsets/ER_051023_Visvader-plus-PDX_MASTER_CancerEROnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#annot <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Annotation_Porting_091423/050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_pseudoClaudinLow_Annotation.csv"

## to ER Cancer Only - SCSubtype
#outFile <- "050823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_simpleMerge_LogNorm_SCSubtype_ported-to-ERCancerOnly.csv"
#targetlib <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/CancerOnlySubsets/ER_051023_Visvader-plus-PDX_MASTER_CancerEROnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#annot <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/SCSubtyper_testing/062823_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_SCSubtypes_Annotation.csv"

#### 11/6/23 Runs
#setwd("~/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/CancerOnlySubsets/Oct2023_Merges")
#sourcelib <- "101023_rerun_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#outDir <- "./"

## to ER Cancer Only - SCSubtype
#outFile <- "ER_101023_Visvader-plus-PDX_MASTER_CancerEROnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_SCSubtypes_Annotation_PortedFrom-101023CancerOnly.csv"
#targetlib <- "ER_101023_Visvader-plus-PDX_MASTER_CancerEROnly_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_LibraryID.csv"
#annot <- "101023_rerun_Visvader-plus-PDX_MASTER_CancerOnly_SimpleMerge_PreLogNorm_SCSubtypes_Annotation.csv"

sourcelib <- opt$sourcelib
outDir <- paste0(opt$outdir,"/")
outFile <- paste0(opt$outdir,"/", opt$outfile)
targetlib <- opt$targetlib
annot <- opt$annotation

print("Summary of input options:\n")
print(paste("Source Library: ", sourcelib))
print(paste("Target Library: ", targetlib))
print(paste("Output Directory and File:" , outFile))
print(paste("Annotation:" , annot))

print("Porting annotations...")
setwd(outDir)

lib_s <- read.csv(sourcelib, header=TRUE)
lib_t <- read.csv(targetlib, header=TRUE)
annot_s <- read.csv(annot, header=TRUE)

## Extract numbers from barcodes
f2 <- function(x,n){return(unlist(strsplit(x,split='-'))[n])}

lib_s$code <- unlist(lapply(lib_s$barcode, f2, n=1))
lib_s$num <- unlist(lapply(lib_s$barcode, f2, n=2))

lib_t$code <- unlist(lapply(lib_t$barcode, f2, n=1))
lib_t$num <- unlist(lapply(lib_t$barcode, f2, n=2))

annot_s$code <- unlist(lapply(annot_s$barcode, f2, n=1))
annot_s$num <- unlist(lapply(annot_s$barcode, f2, n=2))
annot_s_name <- names(annot_s)[2]

## Make a mapping table

map_s <- unique(lib_s[,c("LibraryID","num")])
map_t <- unique(lib_t[,c("LibraryID","num")])

map <- merge(map_s, map_t, all.x = FALSE, all.y = TRUE, by = "LibraryID", sort=FALSE)

map_a <- merge(annot_s, map, all = FALSE, by.x = "num", by.y = "num.x", sort=FALSE)

ported <- data.frame(barcode = paste(map_a$code,map_a$num.y,sep="-"), annotation = map_a[,annot_s_name])
names(ported) <- c("barcode", paste0("ported_",annot_s_name))

## Save results
write.csv(ported, file=outFile, row.names = FALSE, quote = FALSE)








