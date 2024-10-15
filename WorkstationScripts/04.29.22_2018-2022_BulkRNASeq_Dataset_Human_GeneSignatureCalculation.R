## Amy Olex
## 04/29/22
## Script to calculate gene signature scores based off of lists of genes

library("AnnotationDbi")
library(dplyr)

setwd("~/Google Drive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq")

## load in data matrix
data <- read.delim("2018-2022_merged_data_Human/2022.04.29_Log2TPM_Merged_BulkRNASeq_199Samples_noPt_noMouseNormal_noNormalization_DATA.txt", header=TRUE)
samples <- read.delim("2018-2022_merged_data_Human/2022.04.29_Log2TPM_Merged_BulkRNASeq_199Samples_noPt_noMouseNormal_noNormalization_METADATA.txt", header=TRUE)

ensemble_ids <- row.names(data)

#load in gene signature file
sigs <- read.csv("869_gene_signatures_Entrez_Gene_ID.csv", header=TRUE)

gene_sig_mtx <- matrix(0, nrow = length(names(sigs)), ncol = length(names(data)))
colnames(gene_sig_mtx) <- names(data)
rownames(gene_sig_mtx) <- names(sigs)

for(s in rownames(gene_sig_mtx)){
  #get first list of entrez IDs
  entrez_ids <- sigs[!is.na(sigs[,s]),s]
  #identify Ensemble IDs (ususally a many to 1 mapping, so this list will be longer)
  ensembl<-unlist(mapIds(org.Hs.eg.db, keys=as.character(entrez_ids), column='ENSEMBL', keytype='ENTREZID', multiVals='list'))
  #extract data rows that are in the list
  sig_data <- data[row.names(data) %in% ensembl,,drop=FALSE]
  #per sample signature score
  x <- colMeans(sig_data)
  #add to matrix
  gene_sig_mtx[s,] <- x
}

gene_sig_mtx_wTotal <- as.data.frame(gene_sig_mtx)
gene_sig_mtx_wTotal$TotalSigScore <- rowMeans(gene_sig_mtx_wTotal)
gene_sig_mtx_wTotal$SignatureName <- row.names(gene_sig_mtx_wTotal)

gene_sig_mtx_wTotal <- relocate(gene_sig_mtx_wTotal, c("SignatureName","TotalSigScore"))


write.table(gene_sig_mtx_wTotal, file="2022.04.29_Harrell_Merged_Human_869_GeneSignatureScores.tsv", row.names = FALSE, quote = FALSE, sep="\t")


