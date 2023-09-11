## Psudobulk PAM50 subtyping for scRNASeq data
## Amy Olex 
## 9/11/23

library("optparse")
library("Seurat")
library("genefu")
library("DESeq2")
library("dplyr")
library("org.Hs.eg.db")
library("AnnotationDbi")
#library("xtable")
#library("rmeta")
#library("Biobase")
#library("caret")

############
## Helper functions
###########

getGenefuEntrezAnnots <- function(my_keys){
  ### get Entrez IDs
  gene.list <- mapIds(org.Hs.eg.db, keys = my_keys, column = "ENTREZID", keytype = "SYMBOL")
  ## remove NA's
  gene.list.filt <- gene.list[which(!is.na(gene.list))]
  ## convert to dataframe and add columns compatable with genefu
  gene.list.filt <- as.data.frame(gene.list.filt)
  names(gene.list.filt) <- "EntrezGene.ID"
  gene.list.filt$Gene.Symbol <- row.names(gene.list.filt)
  
  return(gene.list.filt)
}

formatPAM50Predictions <- function(preds){
  
  preds[preds == 0] <- NA
  melt_preds <- melt(preds, na.rm=TRUE)
  row.names(melt_preds) <- melt_preds$Var1
  melt_preds <- melt_preds[,"Var2", drop=FALSE]
  names(melt_preds) <- "PAM50"
  
  return(melt_preds)
}


options(future.globals.maxSize = 100000 * 1024^2)

start_time <- Sys.time()

#infile = "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Pseudobulk_PAM50_subtyping/100422_Visvader-MiniSubset_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
#WD = "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/Harrell_LocalTesting/Pseudobulk_PAM50_subtyping/"
#setwd(WD)

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-i", "--infile"), type="character", 
              help="Path and file name of the input RData file with the merged sc samples.", metavar="character"),
  make_option(c("-g", "--groupby"), type="character", 
              help="Seurat identity you want to group by to generate pseudobulk data. Default is orig.ident.",
              default = "orig.ident", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("An input RData file with Seurat single cell object must be provided.", call.=FALSE)
}


runID <- opt$runid
outDir <- paste0(opt$outdir,"/")
inFile <- opt$infile
groups <- opt$groupby

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("RData File: ", inFile))
print(paste("Output Directory:" , outDir))
print(paste("Group By:" , groups))

print("Starting Pseudo PAM50 Subtype Analysis...")
setwd(outDir)

#inFile <- "100422_Visvader-MiniSubset_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RDS"

## Import Seurat Object
seurat.obj <- load(inFile)

## Generate pseudobulk gene expression data
pseudobulk <- AggregateExpression(seurat.obj, assays="RNA", group.by = groups)

## create metadata dataframe
meta_data <- data.frame(sampleID = colnames(pseudobulk$RNA), row.names = colnames(pseudobulk$RNA))

## Normalize the read count data using DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(pseudobulk$RNA, digits = 0), colData = meta_data, design = ~ 1)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
#dds <- estimateDispersions(dds)
#vst_counts <- getVarianceStabilizedData(dds)
#head(vst_counts)

## Convert Gene Symbols to Entrez IDs and generate a gene annotation matrix. This also remove any gene symbol that does not have an ENtrexID in the database.

gene.list.filt <- getGenefuEntrezAnnots(my_keys = row.names(normalized_counts))

## filter normalized counts to only those rows that are still in the filtered gene list.
normalized_counts_filt <- normalized_counts[row.names(gene.list.filt),]

## sanity check when debugging
## There are 47 of the 50 genes for the PAM50 signature availiable in the data set:
#length((which(rownames(normalized_counts_filt) %in% rownames(pam50.robust$centroids.map))))


## run PAM50 classification
pam50_predictions_all <- molecular.subtyping(sbt.model = 'pam50', data = t(normalized_counts_filt), annot = gene.list.filt, do.mapping = FALSE)

pam50_predictions_all$subtype.crisp
pam50_predictions_all$subtype.proba

predictions <- formatPAM50Predictions(preds = pam50_predictions_all$subtype.crisp)

## Merge predictions with metadata
meta_data <- merge(meta_data, predictions, by="row.names")
row.names(meta_data) <- meta_data$Row.names
meta_data <- meta_data[,-1]



## run ClaudinLow classification (https://rdrr.io/bioc/genefu/man/claudinLow.html)
## Using code from https://github.com/clfougner/ClaudinLow/blob/master/Code/METABRIC_patientData.r#L111

## import claudin data
data(claudinLowData)

print("Identifying claudin-low tumors using the nine-cell line predictor (Prat et al. 2010)")

# Find entrez IDs for claudin low genes and identify those available in input scData data
entrezID_CLgenes <- claudinLowData$fnames
overlappingCL_entrezID <- gene.list.filt[which(gene.list.filt$EntrezGene.ID %in% entrezID_CLgenes),]

# Select relevant rows - I don't think I need this because I filter on line 170
#normalized_counts_filtCL <- normalized_counts_filt[row.names(normalized_counts_filt) %in% overlappingCL_entrezID$Gene.Symbol, ]

# Train centroids based on available genes
trainingData <- claudinLowData
trainingData$xd <- medianCtr(trainingData$xd)
trainingData$xd <- trainingData$xd[rownames(trainingData$xd) %in% overlappingCL_entrezID$EntrezGene.ID, ]

# Scale normalized data and training data
normalized_counts_filt_scaled <- t(scale(t(normalized_counts_filt), scale = TRUE, center = TRUE))
normalized_counts_filt_scaled <- normalized_counts_filt_scaled[rownames(overlappingCL_entrezID), ]
row.names(normalized_counts_filt_scaled) <- overlappingCL_entrezID$EntrezGene.ID

trainingData$xd <- t(scale(t(trainingData$xd), scale = TRUE, center = TRUE))

# Determine claudin-low status
cl_class <- claudinLow(x = trainingData$xd, classes = as.matrix(trainingData$classes$Group, ncol = 1), y = normalized_counts_filt_scaled, distm = "euclidean")
predsCL <- cl_class$predictions[,"Call", drop=FALSE]
names(predsCL) <- c("claudinLow")

## Merge predictions with metadata
meta_data <- merge(meta_data, predsCL, by="row.names")
row.names(meta_data) <- meta_data$Row.names
meta_data <- meta_data[,-1]


###### Now add in the annotations for the PAM50 subtyping back into the Seurat object as well as printing out the annotation CSV
print("Adding annotations to Seurat Obj...")

new.annots <- as.data.frame(seurat.obj$orig.ident)
new.annots$barcode <- row.names(new.annots)
names(new.annots) <- c("id", "barcodes")

new.annots.merged <- merge(new.annots, meta_data, by.x = "id", by.y = "row.names", sort = FALSE)
row.names(new.annots.merged) <- new.annots.merged$barcodes

pam50.annots <- new.annots.merged[,c("barcodes","PAM50"), drop=FALSE]
CL.annots <- new.annots.merged[,c("barcodes","claudinLow"), drop=FALSE]

seurat.obj$PAM50 <- new.annots.merged[,"PAM50"]
seurat.obj$claudinLow <- new.annots.merged[,"claudinLow"]

print("Saving data...")
##Writing out output files (rownames remain the same for both)
write.table(pam50.annots, paste0(outDir,runID,"_pseudoPAM50_Annotation.csv"), sep=",", quote = FALSE, row.names = FALSE)
write.table(CL.annots, paste0(outDir,runID,"_pseudoClaudinLow_Annotation.csv"), sep=",", quote = FALSE, row.names = FALSE)

saveRDS(seurat.obj, file = paste0(outDir,runID,"_PseudoPAM50_seuratObj.RDS"))

print("Pseudo PAM50 subtype analysis completed!")












########## The following was code I was playing with, but the method is not robust to implement generalizably, so will save it for another time.





### Using the Core Claudin Low gene signature from the following GitHub:
### https://github.com/clfougner/ClaudinLow/blob/master/Code/METABRIC_CoreClaudinLow.r

#library(ComplexHeatmap)
#library(circlize)
#library(sigclust)
#library(ggsci)

## Read patientData
#exprs <- normalized_counts_filt ##read.table(file = "./Output/METABRIC_patientData.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## Selected claudin-low related genes
#clGenes <- read.table(file = "./ClaudinLowGenes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## Colors
#subtypeColors <- c("#DC362A",
#                   "#ED1E78",
#                   "#333A8E",
#                   "#2F6DB9",
#                   "#267255")

#claudinLowColor <- "yellow"

#uninterestingColor <- "black"


## Subset the selected claudin-low genes
#exprs_CL <- exprs[as.character(clGenes$Hugo_Gene), ]

## Rotate (for scaling)
#exprs_rotated <- t(exprs_CL)

## Center and scale the data
#exprs_centered_scaled <- scale(exprs_rotated, scale = TRUE, center = TRUE)

## Return to correct format for heatmap
#exprs_centered_scaled <- t(exprs_centered_scaled)

### Find core claudin-low samples
#dists <- dist(x = t(exprs_centered_scaled), method = "euclidean")
#clusterObj <- hclust(d = dists, method = "complete")
#whichClust <- cutree(clusterObj, k = 2)
### He randomly selects cluster 2.  Why is it cluster 2??  This is not clear.
### Also he forces it to be in 2 clusters, so even if none of them actually are claudin low they are still forced to be in 2 clusters.
### I don't think I will implement this method right away.
#clClust <- whichClust[whichClust == 2] 
#coreCLsamples <- names(clClust)





