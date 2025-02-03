## Bulk RNASeq PAM50 subtyping
## Amy Olex
## 2/15/2024

library("optparse")
library("genefu")
library("DESeq2")
library("dplyr")
library("org.Hs.eg.db")
library("AnnotationDbi")
data(pam50)
data(pam50.robust)
data(claudinLowData)

############
## Helper functions
###########

getGenefuEntrezAnnots <- function(my_keys){
  ### get Entrez IDs
  gene.list1 <- mapIds(org.Hs.eg.db, keys = my_keys, column = "ENTREZID", keytype = "ENSEMBL")
  gene.list2 <- mapIds(org.Hs.eg.db, keys = my_keys, column = "SYMBOL", keytype = "ENSEMBL")

  gene.list <- merge(as.data.frame(gene.list1), as.data.frame(gene.list2), by="row.names", all=FALSE)

  ## remove NA's
  gene.list.filt <- gene.list[which(!(is.na(gene.list$gene.list1) & is.na(gene.list$gene.list2))),]
  names(gene.list.filt) <- c("Ensembl.ID","EntrezGene.ID", "Gene.Symbol")
  row.names(gene.list.filt) <- gene.list.filt$Ensembl.ID

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

## Function provided by Microsoft Copilot, February 2025.
split_if_needed <- function(x) {
  if (grepl("_", x)) {
    return(strsplit(x, "_")[[1]][1])
  } else {
    return(x)
  }
}


options(future.globals.maxSize = 100000 * 1024^2)

start_time <- Sys.time()


localtest = FALSE
###########################################
#### Local Testing Block
if(localtest){
  setwd("./")
  runID <- "TEST"
  countFile2 <- "/lustre/home/harrell_lab/bulkRNASeq/05_featureCountMatrix/BulkRNASeq_AllGenes_06.20.24_Human_counts.tsv"
  #sampleFile <- "/lustre/home/harrell_lab/bulkRNASeq/config/PAM50_SampleSheet_062024.csv"
  countFile <- "/lustre/home/harrell_lab/bulkRNASeq/05_featureCountMatrix/25.01.31_VCU-BC_PAM50_RSEM_counts.tsv"


  outDir <- "./"
  savedir <- paste0(outDir,runID)
  options(future.globals.maxSize = 3000 * 1024^2)
}
###########################################



option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL,
              help="Required. A unique name for this analysis.", metavar="character"),

  make_option(c("-c", "--countfile"), type="character",
              help="Path and file name of the input raw count data. The first column should be named geneSym and contain the gene symbols for each row.", metavar="character"),

  make_option(c("-s", "--samplefile"), type="character",
              help="Path and file name of the input sample data. Generally retrieved from the master sample sheet. Should at least have the following columns: Sample.ID, Tissue, Treatment, Percent.Human.", metavar="character"),

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

if (is.null(opt$countfile)){
  print_help(opt_parser)
  stop("An input tsv file with raw count values must be provided.", call.=FALSE)
}

if (is.null(opt$samplefile)){
  print_help(opt_parser)
  stop("An input tsv file with sample metadata must be provided.", call.=FALSE)
}


runID <- opt$runid
outDir <- paste0(opt$outdir,"/")
countFile <- opt$countfile
sampleFile <- opt$samplefile


print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("Count File: ", countFile))
print(paste("Sample File: ", sampleFile))
print(paste("Output Directory:" , outDir))


# Load in data
raw_counts <- read.delim(countFile, row.names=1)
formatted_rownames <- sapply(row.names(raw_counts), split_if_needed)
row.names(raw_counts) <- formatted_rownames
names(raw_counts) <- make.names(names(raw_counts), allow_ = FALSE)
genes <- raw_counts[,"geneSym",drop=FALSE]
raw_counts <- raw_counts[,-1]

### Removing SampleSheet Requirement
### I don't need it and with the new non-breast PDXs coming in it will disrupt this workflow.
#samps <- read.delim(sampleFile, row.names=1, sep = ",")
###subset the samples
#samples <- samps[row.names(samps) %in% names(raw_counts),]

#all(row.names(samples)==names(raw_counts))

# create a sample sheet
samples <- data.frame(PAM50 = rep(0, dim(raw_counts)[2]), row.names = names(raw_counts))


# remove rows with zero expression.
raw_counts_nozero <- raw_counts[rowSums(raw_counts) >0,]
genes_nozero <- genes[row.names(raw_counts_nozero),,drop=FALSE]


## Normalize the read count data using DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(raw_counts_nozero, digits = 0), colData = samples, design = ~ 1)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)





## Convert Gene Symbols to Entrez IDs and generate a gene annotation matrix. This also remove any gene symbol that does not have an ENtrezID in the database.

gene.list.filt <- getGenefuEntrezAnnots(my_keys = row.names(normalized_counts))

## filter normalized counts to only those rows that are still in the filtered gene list.
normalized_counts_filt <- normalized_counts[row.names(gene.list.filt),]

## make doubly sure the gene list and normalized counts are in the same order
all(row.names(normalized_counts_filt)==row.names(gene.list.filt))
## change the rownames of the counts matrix to the gene names
row.names(normalized_counts_filt) <- gene.list.filt$Gene.Symbol


print("Running PAM50 Subtyping...")

## run PAM50 classification
pam50_predictions_all <- molecular.subtyping(sbt.model = 'pam50', data = t(normalized_counts_filt), annot = gene.list.filt, do.mapping = FALSE)

#crisp <- as.data.frame(pam50_predictions_all$subtype.crisp)
#pam50_predictions_all$subtype.proba

write.table(file=paste0(outDir,runID,"_PAM50-crisp.tsv"), x=pam50_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
write.table(file=paste0(outDir,runID,"_PAM50-prob.tsv"), x=pam50_predictions_all$subtype.proba, quote=FALSE, sep="\t")

all(row.names(samples)==names(pam50_predictions_all$subtype))

##update annotations:
samples$PAM50 <- pam50_predictions_all$subtype
#samples$PAM50[crisp$Basal==1] <- "Basal"
#samples$PAM50[crisp$Her2==1] <- "Her2"
#samples$PAM50[crisp$LumA==1] <- "LumA"
#samples$PAM50[crisp$LumB==1] <- "LumB"
#samples$PAM50[crisp$Normal==1] <- "Normal"


################ Claudin Low (same method as used for single cell)
print("Running ClaudinLow Subtyping...")
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
#normalized_counts_filt_scaled <- t(scale(t(normalized_counts_filt), scale = TRUE, center = TRUE))
normalized_counts_filt_scaled <- t(medianCtr(t(normalized_counts_filt)))


normalized_counts_filt_scaled <- normalized_counts_filt_scaled[overlappingCL_entrezID$Gene.Symbol,]
row.names(normalized_counts_filt_scaled) <- overlappingCL_entrezID$EntrezGene.ID

#trainingData$xd <- t(scale(t(trainingData$xd), scale = TRUE, center = TRUE))

# Determine claudin-low status
cl_class <- claudinLow(x = trainingData$xd, classes = as.matrix(trainingData$classes$Group, ncol = 1), y = normalized_counts_filt_scaled, distm = "euclidean")
predsCL <- cl_class$predictions[,"Call", drop=FALSE]
names(predsCL) <- c("claudinLow")

## Merge predictions with metadata
all(row.names(samples)==row.names(predsCL))
samples$ClaudinLow <- predsCL$claudinLow
samples$Sample.ID <- row.names(samples)


## save distances
write.table(cl_class$distances, paste0(outDir,runID,"_pseudoClaudinLow_Distances.tsv"), sep="\t", quote = FALSE, row.names = TRUE)


write.table(file=paste0(outDir,runID,"_SubtypedSamples.tsv"), x=samples, quote=FALSE, sep="\t", row.names = FALSE)

print("COMPLETED!")






################## Old ClaudinLow Method
#print("Running ClaudinLow Subtyping...")
## run ClaudinLow classification
#clow_predictions_all <- molecular.subtyping(sbt.model = 'claudinLow', data = t(normalized_counts_filt), annot = gene.list.filt, do.mapping = TRUE)

#crisp_clow <- as.data.frame(clow_predictions_all$subtype.crisp)
#pam50_predictions_all$subtype.proba

#write.table(file=paste0(runID,"_ClaudinLow-crisp.tsv"), x=clow_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
#write.table(file=paste0(runID,"_ClaudinLow-prob.tsv"), x=clow_predictions_all$subtype.proba, quote=FALSE, sep="\t")

#all(row.names(samples)==row.names(crisp_clow))

##update annotations:
#samples$ClaudinLow[crisp_clow$Claudin==1] <- "ClaudinLow"
#samples$ClaudinLow[crisp_clow$Others==1] <- "Others"

#samples$Sample.ID <- row.names(samples)
## write out updated sample file






