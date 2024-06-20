## Amy Olex
## 5/28/24
## HiCAT Immune CEll Typing dara conversion

library(reticulate)
library(Seurat)
library(sceasy)
library(anndata)
use_condaenv('env_py3.8')
loompy <- reticulate::import('loompy')


##converting the Seurat data to anndata

#load("~/Desktop/CCTR_LOCAL_Analysis_noBackups/Harada-Deb_P01_LocalOnly/Deb_Project2/SimpleMerges/SimpleMerge_Mouse_noForce_Log2_GRCh38_052424_Seurat_simpleMerge_LogNormalize_Annotated.RData")
#outFile='SimpleMerge_Mouse_noForce_Log2_GRCh38_052424_Seurat_simpleMerge_LogNormalize_Annotated.h5ad'


setwd("~/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/HiCAT_CellTyping")
inFile <- readRDS("/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/SignatureScoring_Testing/PyMT_Primary_Master_merge_withLuc/031924_PrimaryTumor_pooled_PyMT_Master_SimpleMerge_mm10firefly_ambiantadj_withLuc_Seurat_simpleMerge_LogNormalize_Annotated.rds")
outFile='031924_PrimaryTumor_pooled_PyMT_Master_SimpleMerge_mm10firefly_ambiantadj_withLuc_Seurat_simpleMerge_LogNormalize_Annotated.h5ad'


# seurat to anndata
sceasy::convertFormat(inFile, from="seurat", to="anndata",
                      outFile=outFile)

# anndata to seurat
#sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
#                      outFile='filename.rds')


# Running HiCAT in R
### Load MarkerCount python wrapper
mkrcnt <- import("MarkerCount.hicat")

## get data
adata_test <- read_h5ad(outFile)
X_test <- adata_test$to_df()
mkr_file <- '/Users/alolex/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/SingleCell/HiCAT_cell_markers_rndsystems_mm.tsv'

## run HiCAT
lst_res <- mkrcnt$HiCAT( X_test, marker_file = mkr_file, log_transformed = FALSE ) 

## get results
df_pred <- lst_res[[1]]
df_pred$barcode <- row.names(df_pred)
summary <- lst_res[[2]]

## save results
write.csv(df_pred[,c("barcode","cell_type_major")], file="031924_PrimaryTumor_pooled_PyMT_Master_SimpleMerge_mm10firefly_ambiantadj_withLuc_Seurat_simpleMerge_LogNormalize_MajorCellType.csv", row.names=FALSE, quote = FALSE)
write.csv(df_pred[,c("barcode","cell_type_minor")], file="031924_PrimaryTumor_pooled_PyMT_Master_SimpleMerge_mm10firefly_ambiantadj_withLuc_Seurat_simpleMerge_LogNormalize_MinorCellType.csv", row.names=FALSE, quote = FALSE)
write.csv(df_pred[,c("barcode","cell_type_subset")], file="031924_PrimaryTumor_pooled_PyMT_Master_SimpleMerge_mm10firefly_ambiantadj_withLuc_Seurat_simpleMerge_LogNormalize_MinorSubsetCellType.csv", row.names=FALSE, quote = FALSE)
write.csv(df_pred, file="031924_PrimaryTumor_pooled_PyMT_Master_SimpleMerge_mm10firefly_ambiantadj_withLuc_Seurat_simpleMerge_LogNormalize_HiCAT-Prediction.csv", row.names=FALSE, quote = FALSE)






