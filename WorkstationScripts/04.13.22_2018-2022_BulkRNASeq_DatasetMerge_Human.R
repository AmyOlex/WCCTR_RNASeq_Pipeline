## Amy Olex
## 04/13/22
## Script to merge the bulk RNASeq data from 2018 and 2022 into a single dataset

setwd("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq/2018-2022_merged_data_Human")

#require("hgu219.db")
#library(corrplot)
#library(NMF)
#library(RNASeqBits)
library(ggplot2)
#library(genefilter)
library(ggrepel)
#library(WGCNA)
library(sva)


###### Code run on Fenn:
##load data

TPM2022data <- read.delim("TPM_human_AllGenes_2022_BulkRNASeq_04.13.22.tsv", stringsAsFactors = FALSE, header=TRUE, row.names=1)
TPM2021data <- read.delim("TPM_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv", stringsAsFactors = FALSE, header=TRUE, row.names=1)
TPM2017data <- read.delim("TPM_human_AllGenes_DataFreeze_8.18.17_updated_11.27.17.txt", stringsAsFactors = FALSE, header=TRUE, row.names=1)
TPM_TCGAdata <- read.delim("TCGA-BRCA_STAR-Count_TPM.csv", stringsAsFactors = FALSE, header=TRUE, row.names=1, sep=",")
##remove duplicate PAR_Y ensemble IDs
rownames <- unlist(lapply(strsplit(row.names(TPM_TCGAdata), split=".", fixed=TRUE), function(x) x[1]))
dups <- duplicated(rownames)
parY <- row.names(TPM_TCGAdata)[which(dups)]
tpm_parY <- rowSums(TPM_TCGAdata[parY,2:1226])
TPM_TCGAdata <- TPM_TCGAdata[!(row.names(TPM_TCGAdata) %in% parY),]
#remove version
row.names(TPM_TCGAdata) <- unlist(lapply(strsplit(row.names(TPM_TCGAdata), split=".", fixed=TRUE), function(x) x[1]))

#create the metadata
TCGAsamples <- data.frame(row.names=names(TPM_TCGAdata[2:1227]))
TCGAsamples$PDX.line = "TCGA"
TCGAsamples$Treatment = "TCGA"
TCGAsamples$tissue = "MGT"
TCGAsamples$batch = "0"
TCGAsamples$type = "TCGA"


##################
## Import Data
##################
## Import RNASeq TPM data
TPM2022data <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq/processed_data/TPM_human_AllGenes_2022_BulkRNASeq_04.13.22.tsv", stringsAsFactors = FALSE, header=TRUE, row.names=1)
TPM2021data <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2021_BulkRNASeq/processed_data/TPM_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv", stringsAsFactors = FALSE, header=TRUE, row.names=1)
TPM2017data <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/8-18-17 Raw Count Data Freeze/TPM_human_AllGenes_DataFreeze_8.18.17_updated_11.27.17.txt", stringsAsFactors = FALSE, header=TRUE, row.names=1)


## Import Samples
TPM2017samples <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/8-18-17 Raw Count Data Freeze/script and input files/Sample_Metadata_Freeze_8-18-17_updated_12-15-21.txt", header=TRUE, row.names=1, stringsAsFactors = FALSE)

#TPM2017samples <- read.delim("Sample_Metadata_Freeze_8-18-17_updated_12-15-21.txt", header=TRUE, row.names=1, stringsAsFactors = FALSE)

TPM2017samples$Treatment <- "untreated"
TPM2017samples <- TPM2017samples[,c("PDX.line", "Treatment","tissue","batch", "type")]
## filter 2017 to remove patient samples and normal mouse
TPM2017samples <- TPM2017samples[!(TPM2017samples$PDX.line %in% c("Normal", "Patient 1", "Patient 2", "Patient 3")),]
TPM2017data <- TPM2017data[,row.names(TPM2017samples)]


TPM2021samples <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2021_BulkRNASeq/sample_metadata_master_10.28.21.tsv", header=TRUE, row.names=1, stringsAsFactors = FALSE)

#TPM2021samples <- read.delim("sample_metadata_master_10.28.21.tsv", header=TRUE, row.names=1, stringsAsFactors = FALSE)

TPM2021samples$batch = "6"
TPM2021samples$tissue = "MGT"
TPM2021samples$type = "MGT"
TPM2021samples <- TPM2021samples[,c("PDX.line", "Treatment","tissue","batch", "type")]

TPM2022samples <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq/sample_metadata_master_04.13.22.tsv", header=TRUE, row.names=1, stringsAsFactors = FALSE)

#TPM2022samples <- read.delim("sample_metadata_master_04.13.22.tsv", header=TRUE, row.names=1, stringsAsFactors = FALSE)
TPM2022samples <- TPM2022samples[,c("PDX.line", "Treatment","tissue","batch", "type")]


## Merge Data
#TPM <- merge(TPM_TCGAdata, TPM2017data, by="row.names", all=FALSE)
#row.names(TPM) <- TPM$Row.names
#names(TPM) <- make.names(names(TPM))
#TPM <- TPM[,-1]
#TPM <- TPM[,-1]

TPM <- merge(TPM2021data, TPM2017data, by="row.names", all=FALSE)
#TPM <- merge(TPM2021data,TPM, by="row.names", all=FALSE)
#genes <- TPM[,"symbol",drop=FALSE]
row.names(TPM) <- TPM$Row.names
#row.names(genes) <- TPM$Row.names
names(TPM) <- make.names(names(TPM))
TPM <- TPM[,-1]
TPM <- TPM[,-1]


TPM <- merge(TPM2022data, TPM, by="row.names", all=FALSE)
genes <- TPM[,"symbol",drop=FALSE]
row.names(TPM) <- TPM$Row.names
row.names(genes) <- TPM$Row.names
TPM <- TPM[,-1]
TPM <- TPM[,-1]

#samples <- rbind(TPM2022samples,  TPM2021samples, TPM2017samples, TCGAsamples)
samples <- rbind(TPM2022samples,  TPM2021samples, TPM2017samples)

row.names(samples) <- make.names(row.names(samples))
TPMdata <- TPM[,row.names(samples)]

## save data
#write.table(TPMdata, file="Harrell-TCGA_merged_2017-2022_TPM_noNorm_042222.csv", quote=FALSE, row.names=TRUE, sep=",")
#write.table(samples, file="Harrell-TCGA_merged_2017-2022_MetaData_042222.csv", quote=FALSE, row.names=TRUE, sep=",")

                          
################
## Filter low expression from RNASeq and array
################

#TPM <- read.delim("Harrell-TCGA_merged_2017-2022_TPM_noNorm_042222.csv", row.names=1, sep=",")
#samples <- read.delim("Harrell-TCGA_merged_2017-2022_MetaData_042222.csv", row.names=1, sep=",")


## First remove zero sum rows from TPM data.
#TPMdata_nozero <- TPMdata[rowSums(TPMdata[,2:ncol(TPMdata)]) > 0,]
TPMdata_nozero <- TPMdata[rowSums(TPMdata) > 0,]

TPMdata_nozero_log2 <- log2(TPMdata_nozero+1)
## Identify RNASeq 5% quartile
#rnaseq_quant_cutoff <- quantile(as.matrix(TPMdata_nozero), probs=0.05)[1]
## Remove all genes with more than 90% of samples having a value below the cutoff.



#TPMdata_filtered <- TPMdata_nozero[which(rowSums(TPMdata_nozero <= rnaseq_quant_cutoff) < floor(.9*ncol(TPMdata_nozero))),]

#TPMdata_log2 <- log2(TPMdata_filtered+1)

## save log2 data
write.table(TPMdata_nozero_log2, file="Harrell-TCGA_merged_2017-2022_Log2TPM_noZero_noNorm_042222.csv", quote=FALSE, row.names=TRUE, sep=",")
##### I got the above saved.  will need to import, run PCA and create the figures, then normalize and run PAM50 subtyping.

TPMdata_nozero_log2 <- read.delim("Harrell-TCGA_merged_2017-2022_Log2TPM_noZero_noNorm_042222.csv", row.names=1, sep=",")
samples <- read.delim("Harrell-TCGA_merged_2017-2022_MetaData_042222.csv", row.names=1, sep=",")

rowStd <- apply(TPMdata_nozero_log2,1,sd)

top2000_log2 <- TPMdata_nozero_log2[order(rowStd, decreasing=TRUE),][1:2000,]
top5000_log2 <- TPMdata_nozero_log2[order(rowStd, decreasing=TRUE),][1:5000,]
#top2000_log2 <- log2(top2000+1)
###############
## Write out TPM Log2 un-normalized data
###############

## Print out TPM filt file and samples
write.table(file="2022.04.13_Log2TPM_Merged_BulkRNASeq_202Samples_noPt_noMouseNormal_noNormalization_DATA.txt", x=TPMdata_nozero_log2, quote=FALSE, sep="\t", row.names=TRUE)
write.table(file="2022.04.13_Log2TPM_Merged_BulkRNASeq_202Samples_noPt_noMouseNormal_noNormalization_METADATA.txt", x=samples, quote=FALSE, sep="\t", row.names=TRUE)

## remove adjacent normal
samples <- samples[samples$type != "Normal.Adj",]
TPMdata_nozero_log2 <- TPMdata_nozero_log2[,row.names(samples)]

## Print out TPM filt file and samples
write.table(file="2022.04.29_Log2TPM_Merged_BulkRNASeq_199Samples_noPt_noMouseNormal_noNormalization_DATA.txt", x=TPMdata_nozero_log2, quote=FALSE, sep="\t", row.names=TRUE)
write.table(file="2022.04.29_Log2TPM_Merged_BulkRNASeq_199Samples_noPt_noMouseNormal_noNormalization_METADATA.txt", x=samples, quote=FALSE, sep="\t", row.names=TRUE)


## remove genes with low variance.  use low expression instead, below 5 TPM 90% or more
#TPMdata_filtered2 <- as.data.frame(varFilter(as.matrix(TPMdata_selected[,2:ncol(TPMdata_selected)])))

## PCA Analysis before normalization
pca_before <-  prcomp(t(top5000_log2))

to_plot <- data.frame(pca_before$x[,c("PC1","PC2")], batch = samples$batch, pdx = samples$PDX.line, type=samples$type, Sample.ID = row.names(samples), stringsAsFactors = F)
       
#ggplot(data = to_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2)), color = batch) + 
#  geom_point(aes(color = batch), size = 1) + theme_bw() + ggtitle("All samples - colored by batch (batch 6 is the 2021 set)")

png(filename = "PCAplot_Top5000_log2TPM_human_ZeroExpRemoved_NoNorm_byBatch_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "batch"
ggplot(data = to_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()

#ggplot(data = to_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2)), color = pdx) + 
#  geom_point(aes(color = pdx), size = 1) + theme_bw() + ggtitle("All samples - colored by PDX")

png(filename = "PCAplot_Top5000_log2TPM_human_ZeroExpRemoved_NoNorm_byPDX_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "pdx"
ggplot(data = to_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()

#ggplot(data = to_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2)), color = type) + 
#  geom_point(aes(color = type), size = 1) + theme_bw() + ggtitle("All samples - colored by tissue type")


png(filename = "PCAplot_Top5000_log2TPM_human_ZeroExpRemoved_NoNorm_byType_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "type"
ggplot(data = to_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()


##run PAM50 analysis
library(genefu)

tpm_t <- t(TPMdata_nozero_log2)
colnames(tpm_t) <- genes$symbol

## run PAM50 classification
pam50_predictions_all <- molecular.subtyping(sbt.model = "pam50", data = tpm_t, annot = genes, do.mapping = FALSE)

crisp <- as.data.frame(pam50_predictions_all$subtype.crisp)
pam50_predictions_all$subtype.proba

write.table(file="PAM50-crisp_TCGA-Harrell_merge_042222.txt", x=pam50_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
write.table(file="PAM50-prob_TCGA-Harrell_merge_042222.txt", x=pam50_predictions_all$subtype.proba, quote=FALSE, sep="\t")

##update annotations:
samples$PAM50[crisp$Basal==1] <- "Basal"
samples$PAM50[crisp$Her2==1] <- "Her2"
samples$PAM50[crisp$LumA==1] <- "LumA"
samples$PAM50[crisp$LumB==1] <- "LumB"
samples$PAM50[crisp$Normal==1] <- "Normal"

write.table(file="Harrell-TCGA_merged_2017-2022_MetaData_wPAM50_042222.csv", x=samples, row.names=TRUE, quote=FALSE, sep=",")


### Also running on Combat cleaned data to compare:
tpm_t <- t(combat_edata)
colnames(tpm_t) <- genes$symbol

## run PAM50 classification using ComBatNorm data
pam50_predictions_all <- molecular.subtyping(sbt.model = "pam50", data = tpm_t, annot = genes, do.mapping = FALSE)

crisp <- as.data.frame(pam50_predictions_all$subtype.crisp)
pam50_predictions_all$subtype.proba

write.table(file="PAM50-crisp_TCGA-Harrell_merge_ComBat_042222.txt", x=pam50_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
write.table(file="PAM50-prob_TCGA-Harrell_merge_ComBat_042222.txt", x=pam50_predictions_all$subtype.proba, quote=FALSE, sep="\t")

##update annotations:
samples$PAM50_ComBat[crisp$Basal==1] <- "Basal"
samples$PAM50_ComBat[crisp$Her2==1] <- "Her2"
samples$PAM50_ComBat[crisp$LumA==1] <- "LumA"
samples$PAM50_ComBat[crisp$LumB==1] <- "LumB"
samples$PAM50_ComBat[crisp$Normal==1] <- "Normal"

write.table(file="Harrell-TCGA_merged_2017-2022_MetaData_wPAM50_042222.csv", x=samples, row.names=TRUE, quote=FALSE, sep=",")


## run PAM50 classification using only Harrell Samples
tpm_t <- t(TPMdata_nozero_log2[,samples$PDX.line != "TCGA"])
colnames(tpm_t) <- genes$symbol

pam50_predictions_all <- molecular.subtyping(sbt.model = "pam50", data = tpm_t, annot = genes, do.mapping = FALSE)

crisp <- as.data.frame(pam50_predictions_all$subtype.crisp)
pam50_predictions_all$subtype.proba

write.table(file="PAM50-crisp_HarrellOnly_merge_NoNorm_042222.txt", x=pam50_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
write.table(file="PAM50-prob_HarrellOnly_merge_NoNorm_042222.txt", x=pam50_predictions_all$subtype.proba, quote=FALSE, sep="\t")

m <- data.frame(matrix(0, ncol = 5, nrow = 1226))
names(m) <- c("Basal","Her2","LumA","LumB","Normal")
crisp_extended <- rbind(crisp, m)

##update annotations:
samples$PAM50_Harrell[crisp_extended$Basal==1] <- "Basal"
samples$PAM50_Harrell[crisp_extended$Her2==1] <- "Her2"
samples$PAM50_Harrell[crisp_extended$LumA==1] <- "LumA"
samples$PAM50_Harrell[crisp_extended$LumB==1] <- "LumB"
samples$PAM50_Harrell[crisp_extended$Normal==1] <- "Normal"

write.table(file="Harrell-TCGA_merged_2017-2022_MetaData_wPAM50_042222.csv", x=samples, row.names=TRUE, quote=FALSE, sep=",")








to_plot <- data.frame(pca_before$x[,c("PC1","PC2")], batch = samples$batch, pdx = samples$PDX.line, type=samples$type, Sample.ID = row.names(samples), pam50=samples$PAM50, stringsAsFactors = F)


png(filename = "PCAplot_Top5000_log2TPM_human_ZeroExpRemoved_NoNorm_byPAM50_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "pam50"
ggplot(data = to_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()











#samples_filt_pam50 <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2021_BulkRNASeq/2018-2021_merged_data/Annotations_PAM50_151Samples.txt")
#row.names(samples_filt_pam50) <- samples_filt_pam50$Row.names
#all(row.names(samples_filt_pam50) == row.names(samples_filt))

#to_plot_filt2 <- data.frame(pca_before_filt$x[,c("PC1","PC2")], batch = samples_filt_pam50$batch, pdx = samples_filt_pam50$PDX.line, type=samples_filt_pam50$type, PAM50=samples_filt_pam50$PAM50, Sample.ID = row.names(samples_filt_pam50), stringsAsFactors = F)

#png(filename = "PCAplot_log2TPM_human_ZeroExpRemoved_NoNorm_byPAM50Subtype_FILTERED_12.15.21.png", res = 150, width=2000, height=1000)  

#colorby <- "PAM50"
#ggplot(data = to_plot_filt2, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
#  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
#  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
#  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
#  geom_text_repel(colour = "black", size = 3) +
#  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
#  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
#  labs(color = colorby) +
#  theme_bw() +
#  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
#  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

#dev.off()

############
## Batch Effect Correction
############

modcombat <- model.matrix(~1, data = samples[,"batch",drop=FALSE])
combat_edata <- ComBat(dat = TPMdata_nozero_log2, batch = samples$batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
combat_seqdata <- ComBat_seq(counts = TPMdata_nozero_log2, batch = samples$batch)

#combat_towrite <- merge(genes, combat_edata, by="row.names")
write.table(file="ComBatSeqNorm_Harrell-TCGA_merged_2017-2022_Log2TPM_noZero_noNorm_042222.csv", x=combat_edata, quote=FALSE, sep=",", row.names=TRUE)

#rowStd <- apply(combat_edata,1,sd)
top2000_norm <- combat_edata[order(rowStd, decreasing=TRUE),][1:2000,]
top5000_norm <- combat_seqdata[order(rowStd, decreasing=TRUE),][1:5000,]


## PCA Analysis of ComBat normalized merged data
pca_combat <-  prcomp(t(top5000_norm))

to_plot_combat <- data.frame(pca_combat$x[,c("PC1","PC2")], batch = samples$batch, pdx = samples$PDX.line, type=samples$type, pam50=samples$PAM50, Sample.ID = row.names(samples), stringsAsFactors = F)

png(filename = "PCAplot_top5000_log2TPM_human_ZeroExpRemoved_ComBatSeqNorm_byBatch_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "batch"
ggplot(data = to_plot_combat, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_combat)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_combat)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()

png(filename = "PCAplot_top5000_log2TPM_human_ZeroExpRemoved_ComBatSeqNorm_byPDX_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "pdx"
ggplot(data = to_plot_combat, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_combat)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_combat)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()

png(filename = "PCAplot_top5000_log2TPM_human_ZeroExpRemoved_ComBatSeqNorm_byType_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "type"
ggplot(data = to_plot_combat, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_combat)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_combat)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()


png(filename = "PCAplot_top5000_log2TPM_human_ZeroExpRemoved_ComBatSeqNorm_byPAM50_04.22.22.png", res = 150, width=2000, height=1000)  

colorby <- "pam50"
ggplot(data = to_plot_combat, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_combat)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_combat)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()

#to_plot_combat2 <- data.frame(pca_combat$x[,c("PC1","PC2")], batch = samples_filt_pam50$batch, pdx = samples_filt_pam50$PDX.line, type=samples_filt_pam50$type, PAM50=samples_filt_pam50$PAM50, Sample.ID = row.names(samples_filt_pam50), stringsAsFactors = F)

#png(filename = "PCAplot_log2TPM_human_ZeroExpRemoved_ComBatNorm_byPAM50Subtype_FILTERED_12.15.21.png", res = 150, width=2000, height=1000)  

#colorby <- "PAM50"
#ggplot(data = to_plot_combat2, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
#  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
#  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
#  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
#  geom_text_repel(colour = "black", size = 3) +
#  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
#  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
#  labs(color = colorby) +
#  theme_bw() +
#  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
#  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

#dev.off()

########################## DID NOT use code below this line for the 2022 merge ################################
############
## Extracting UCD52 Only
############


UCD52_log2 <- TPMdata_log2[,samples$PDX.line%in%c("UCD52","PT52")]
UCD52_samples <- samples[samples$PDX.line%in%c("UCD52","PT52"),]

pca_ucd52 <-  prcomp(t(UCD52_log2))
ucd52_plot <- data.frame(pca_ucd52$x[,c("PC1","PC2")], batch = UCD52_samples$batch, pdx = UCD52_samples$PDX.line, type=UCD52_samples$type, sample.id=row.names(UCD52_samples), stringsAsFactors = F)

#ggplot(data = ucd52_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2)), color = batch) +
#  geom_point(aes(color = batch), size = 1) + theme_bw()

png(filename = "PCAplot_UCD52_log2TPM_human_ZeroExpRemoved_NoNorm_byBatchWID_12.15.21.png", res = 150, width=1500, height=600)  

colorby <- "batch"
ggplot(data = ucd52_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = sample.id)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()

png(filename = "PCAplot_UCD52_log2TPM_human_ZeroExpRemoved_NoNorm_byBatch_12.15.21.png", res = 150, width=800, height=600)  

colorby <- "batch"
ggplot(data = ucd52_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()




png(filename = "PCAplot_UCD52_log2TPM_human_ZeroExpRemoved_NoNorm_byType_12.15.21.png", res = 150, width=700, height=600)  

colorby <- "type"
ggplot(data = ucd52_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2), label = pdx)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_before)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_before)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()






ggplot(data = ucd52_plot, aes(x = as.numeric(PC1), y = as.numeric(PC2)), color = type) +
  geom_point(aes(color = type), size = 1) + theme_bw()
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  geom_hline(yintercept = 0, colour = "darkgrey", size = .25) +
  geom_vline(xintercept = 0, colour = "darkgrey", size = .25) +
  ggtitle(paste("PCA with batch, coloring by ", colorby)) +
  
  scale_color_manual(values=c("red", "pink", "burlywood3", "darkblue","gold", "cyan", "darkgreen","green","lightgreen", "orange", "magenta","purple","blue","black","lightgrey")) +
  #geom_text_repel(colour = "black", size = 3, point.padding = unit(0.25, "lines")) +
  labs(color = colorby) +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca)$importance[2,2] * 100, digits = 2), "% variability" )) +
  

png(filename=paste("PCA_combat_Breast.png", sep=""), width=2000, height=1200, res=250)
plot(pt)
dev.off()








## COMBAT normalization
modcombat <- model.matrix(~1, data = common_annot[,"batch",drop=FALSE])
combat_edata <- ComBat(dat = combined_dat, batch = common_annot$batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

combat_towrite <- merge(common_geneIDs, combat_edata, by="row.names")
write.table(file="2019.01.31_ComBatNorm_Merged_BreastCancerCell_Microarray-RNASeq.txt", x=combat_towrite, quote=FALSE, sep="\t", row.names=FALSE)


plot(combat_edata[,2],combat_edata[,94])



## PCA Analysis of ComBat normalized merged data
pca <-  prcomp(t(combat_edata))

## The percent variance Explained:
data.frame(summary(pca)$importance)[, 1:min(5, ncol(summary(pca)$importance))]

common_annot_rev <- common_annot[rev(1:nrow(common_annot)),]

colorby <- "type"
pt <- ggplot(data = data.frame(pca$x[rev(1:nrow(pca$x)),], common_annot_rev, samples = common_annot_rev$type, stringsAsFactors = F), 
             aes(x = as.numeric(PC1), y = as.numeric(PC2), label = type)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  geom_hline(yintercept = 0, colour = "darkgrey", size = .25) +
  geom_vline(xintercept = 0, colour = "darkgrey", size = .25) +
  ggtitle(paste("PCA with batch, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 1) +
  scale_color_manual(values=c("red", "pink", "burlywood3", "darkblue","gold", "cyan", "darkgreen","green","lightgreen", "orange", "magenta","purple","blue","black","lightgrey")) +
  #geom_text_repel(colour = "black", size = 3, point.padding = unit(0.25, "lines")) +
  labs(color = colorby) +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca)$importance[2,2] * 100, digits = 2), "% variability" )) +
  theme_bw()

png(filename=paste("PCA_combat_Breast.png", sep=""), width=2000, height=1200, res=250)
plot(pt)
dev.off()





cormat <- cor(combat_edata, method="pearson")
centered_colors <- center.palette(cormat, palette_length = 100)



jpeg(filename=paste("2019.01.31_PDX-BreastCellLine_Correlation.jpeg", sep=""), width=3000, height=3000, res=300)
corrplot(cormat, cl.lim = c(.65, 1), tl.cex=.3, method="color", type="upper", order="hclust", hclust.method="average", col=colorRampPalette(c("blue","blue","blue","blue","blue","white","red"))(200), tl.col="black")
dev.off()

jpeg(filename=paste("2019.01.31_PDX-BreastCellLine_CorrelationHeatmap.jpeg", sep=""), width=6000, height=8000, res=1000)
aheatmap(cormat, distfun="pearson", hclustfun="average", scale="none", color=colorRampPalette(c("white","red"))(200), fontsize=10, cexRow=1, annCol = common_annot$batch)
dev.off()





#### Found plotly code from : https://moderndata.plot.ly/principal-component-analysis-cluster-plotly/
library(plotly)
p_data <- as.data.frame(pca$x)

write.table(file="2019.01.31_PCARevResults_BreastCancerCell_Microarray-RNASeq.txt", x=p_data, quote=FALSE, sep="\t", row.names=TRUE)
write.table(file="2019.01.31_Annots_BreastCancerCell_Microarray-RNASeq.txt", x=common_annot, quote=FALSE, sep="\t", row.names=TRUE)





p <- plot_ly(p_data[94:1110,], type = "scatter", x = as.numeric(p_data[94:1110,"PC1"]) , y = as.numeric(p_data[94:1110,"PC2"]), 
              text = row.names(common_annot)[94:1110],mode = "markers", marker = list(size = 8), color=common_annot$type[94:1110], colors = c("red", "pink", "burlywood3", "darkblue","gold", "cyan", "darkgreen","green","lightgreen", "orange", "magenta","purple","blue","black","lightgrey")) 

p <- add_trace(p, data = p_data[1:93,], x = as.numeric(p_data[1:92,"PC1"]) , y = as.numeric(p_data[1:92,"PC2"]), text = row.names(common_annot)[1:92],
             mode = "markers", marker = list(size = 10), color=common_annot$type[1:92], colors = c("red", "pink", "burlywood3", "darkblue","gold", "cyan", "darkgreen","green","lightgreen", "orange", "magenta","purple","blue","black")) 

p <- layout(p, title = "PCA Plot", 
            xaxis = list(title = "PC 1"),
            yaxis = list(title = "PC 2"))

p








