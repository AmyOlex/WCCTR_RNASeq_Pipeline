## Amy Olex
## 12.15.21
## PAM50 Classification of large Harrell Cohort for MGT samples only.
## 
## UPDATE 11/28/17: Updated to import the latest human percentage adjusted TPM values.
## UPDATE 12/3/17:  Updated to remove the mouse normal tissues from the classification.  Updated to also output the probabilities.
## UPDATE: 5/22/18: Updated to remove all WHIM30 samples.  Also updated to import the TPM values that are not corrected for human percent.
## UPDATE: 7/5/18: Updated to put all WHIM30 samples back in and leave the mouse normal samples out and continue using the TPM values that are not corrected for Human Percent.
## UPDATE: 7/23/18: Updated to remove all samples with less than 1 million human reads.
## UPDATE: 12/15/21: to use new and old merged samples.

library(genefu)

setwd("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2021_BulkRNASeq/2018-2021_merged_data")

## import unadjusted TPM values from data freeze
tpm <- read.delim("2021.12.15_Log2TPM_Merged_BulkRNASeq_153Samples_noPt_noNormal_noNormalization_DATA.txt")

annot <- read.delim("2021.12.15_Log2TPM_Merged_BulkRNASeq_153Samples_noPt_noNormal_noNormalization_METADATA.txt")

gene_annot <- tpm[,c("Row.names","symbol.x")]
names(gene_annot) <- c("Ensemble.ID","Gene.Symbol")
tpm <- tpm[,-1]
tpm <- tpm[,-1]
## Reorder
#tpm <- tpm[,as.character(annot$preferred.sample.name)]
stopifnot(all(as.character(annot$Row.names) == names(tpm)))


############## All Data ##############
## transpose the gene expression matrix
tpm_t <- t(tpm)
colnames(tpm_t) <- gene_annot$Gene.Symbol

## run PAM50 classification
pam50_predictions_all <- molecular.subtyping(sbt.model = "pam50", data = tpm_t, annot = gene_annot, do.mapping = FALSE)

crisp <- as.data.frame(pam50_predictions_all$subtype.crisp)
pam50_predictions_all$subtype.proba

write.table(file="PAM50_151Samples.txt", x=pam50_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
##update annotations:
annot$PAM50[crisp$Basal==1] <- "Basal"
annot$PAM50[crisp$Her2==1] <- "Her2"
annot$PAM50[crisp$LumA==1] <- "LumA"
annot$PAM50[crisp$LumB==1] <- "LumB"
annot$PAM50[crisp$Normal==1] <- "Normal"

write.table(file="Annotations_PAM50_151Samples.txt", x=annot, quote=FALSE, sep="\t")






############## Samples with > 1 million human reads only ##############

## remove all normal samples
annot_filtered_1M <- annot[annot$human.reads >= 1000000,]

## Reorder
tpm_filtered_1M <- tpm[,as.character(annot_filtered_1M$preferred.sample.name)]
stopifnot(all(as.character(annot_filtered_1M$preferred.sample.name) == names(tpm_filtered_1M)))


## transpose the gene expression matrix
tpm_filtered_1M_t <- t(tpm_filtered_1M)
colnames(tpm_filtered_1M_t) <- gene_annot$Gene.Symbol

## run PAM50 classification
pam50_predictions_filtered_1M <- molecular.subtyping(sbt.model = "pam50", data = tpm_filtered_1M_t, annot = gene_annot, do.mapping = FALSE)

pam50_predictions_filtered_1M$subtype.crisp
pam50_predictions_filtered_1M$subtype.proba

############## Samples with > 10 million human reads only ##############

## remove all normal samples
annot_filtered_10M <- annot[annot$human.reads >= 10000000,]

## Reorder
tpm_filtered_10M <- tpm[,as.character(annot_filtered_10M$preferred.sample.name)]
stopifnot(all(as.character(annot_filtered_10M$preferred.sample.name) == names(tpm_filtered_10M)))


## transpose the gene expression matrix
tpm_filtered_10M_t <- t(tpm_filtered_10M)
colnames(tpm_filtered_10M_t) <- gene_annot$Gene.Symbol

## run PAM50 classification
pam50_predictions_filtered_10M <- molecular.subtyping(sbt.model = "pam50", data = tpm_filtered_10M_t, annot = gene_annot, do.mapping = FALSE)

pam50_predictions_filtered_10M$subtype.crisp
pam50_predictions_filtered_10M$subtype.proba







