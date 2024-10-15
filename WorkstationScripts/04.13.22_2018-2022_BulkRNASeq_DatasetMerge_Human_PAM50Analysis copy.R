## Amy Olex
## 04/29/22
## Script to run PAM50 analysis on merged Harrrell data


library(sva)
library(genefu)

setwd("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq/2018-2022_merged_data_Human/")




# Load in data
TPMdata_nozero_log2 <- read.delim("2022.04.29_Log2TPM_Merged_BulkRNASeq_199Samples_noPt_noMouseNormal_noNormalization_DATA.txt", row.names=1)
samples <- read.delim("2022.04.29_Log2TPM_Merged_BulkRNASeq_199Samples_noPt_noMouseNormal_noNormalization_METADATA.txt", row.names=1)
genes <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq/processed_data/Count_human_AllGenes_2022_BulkRNASeq_04.13.22.tsv", header=TRUE, row.names=1, stringsAsFactor=FALSE)
genes <- genes[row.names(TPMdata_nozero_log2),"symbol",drop=FALSE]


# Subset to only MGT samples
samples_mgt <- samples[samples$type == "MGT",]
data_mgt <- TPMdata_nozero_log2[,row.names(samples_mgt)]

##run PAM50 analysis
tpm_t <- t(data_mgt)
colnames(tpm_t) <- genes$symbol

## run PAM50 classification
pam50_predictions_all <- molecular.subtyping(sbt.model = "pam50", data = tpm_t, annot = genes, do.mapping = FALSE)

crisp <- as.data.frame(pam50_predictions_all$subtype.crisp)
pam50_predictions_all$subtype.proba

write.table(file="PAM50-crisp_Harrell-MGTonly_merge_042222.txt", x=pam50_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
write.table(file="PAM50-prob_Harrell-MGTonly_merge_042222.txt", x=pam50_predictions_all$subtype.proba, quote=FALSE, sep="\t")

##update annotations:
samples_mgt$PAM50_MGT[crisp$Basal==1] <- "Basal"
samples_mgt$PAM50_MGT[crisp$Her2==1] <- "Her2"
samples_mgt$PAM50_MGT[crisp$LumA==1] <- "LumA"
samples_mgt$PAM50_MGT[crisp$LumB==1] <- "LumB"
samples_mgt$PAM50_MGT[crisp$Normal==1] <- "Normal"

##read in master annotation file
master <- read.delim("../2017-2022-TCGA_merged_data_Human/Harrell-TCGA_merged_2017-2022_MetaData_wPAM50_042222.csv", row.names=1, header=TRUE, sep=",")
##merge in new PAM50 calls
master_updated <- merge(master,samples_mgt[,"PAM50_MGT",drop=FALSE], all.x=TRUE, by="row.names")
row.names(master_updated) <- master_updated$Row.names
master_updated <- master_updated[,-1]
samples_updated <- merge(samples,samples_mgt[,"PAM50_MGT",drop=FALSE], all.x=TRUE, by="row.names")
row.names(samples_updated) <- samples_updated$Row.names
samples_updated <- samples_updated[,-1]




###run again using MGT and MGTfromMet samples
# Subset to only MGT samples
samples_mgtmet <- samples[samples$type %in% c("MGT","MGTFromMet"),]
data_mgt <- TPMdata_nozero_log2[,row.names(samples_mgtmet)]

##run PAM50 analysis
tpm_t <- t(data_mgt)
colnames(tpm_t) <- genes$symbol

## run PAM50 classification
pam50_predictions_all <- molecular.subtyping(sbt.model = "pam50", data = tpm_t, annot = genes, do.mapping = FALSE)

crisp <- as.data.frame(pam50_predictions_all$subtype.crisp)
pam50_predictions_all$subtype.proba

write.table(file="PAM50-crisp_Harrell-MGT-MGTFromMet-only_merge_042222.txt", x=pam50_predictions_all$subtype.crisp, quote=FALSE, sep="\t")
write.table(file="PAM50-prob_Harrell-MGT-MGTFromMet-only_merge_042222.txt", x=pam50_predictions_all$subtype.proba, quote=FALSE, sep="\t")

##update annotations:
samples_mgtmet$PAM50_MGTFromMet[crisp$Basal==1] <- "Basal"
samples_mgtmet$PAM50_MGTFromMet[crisp$Her2==1] <- "Her2"
samples_mgtmet$PAM50_MGTFromMet[crisp$LumA==1] <- "LumA"
samples_mgtmet$PAM50_MGTFromMet[crisp$LumB==1] <- "LumB"
samples_mgtmet$PAM50_MGTFromMet[crisp$Normal==1] <- "Normal"


##merge in new PAM50 calls
master_updated <- merge(master_updated,samples_mgtmet[,"PAM50_MGTFromMet",drop=FALSE], all.x=TRUE, by="row.names")
samples_updated <- merge(samples_updated,samples_mgtmet[,"PAM50_MGTFromMet",drop=FALSE], all.x=TRUE, by="row.names")




write.table(file="Harrell-TCGA_merged_2017-2022_MetaData_wPAM50_042922.csv", x=master_updated, row.names=FALSE, quote=FALSE, sep=",")
write.table(file="Harrell-merged_2017-2022_MetaData_wPAM50_042922.csv", x=samples_updated, row.names=FALSE, quote=FALSE, sep=",")



