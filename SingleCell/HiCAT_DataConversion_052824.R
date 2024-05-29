## Amy Olex
## 5/28/24
## HiCAT Immune CEll Typing dara conversion

library(reticulate)
library(Seurat)
library(sceasy)

##converting the Seurat data to anndata

load("~/Desktop/CCTR_LOCAL_Analysis_noBackups/Harada-Deb_P01_LocalOnly/Deb_Project2/SimpleMerges/SimpleMerge_Mouse_noForce_Log2_GRCh38_052424_Seurat_simpleMerge_LogNormalize_Annotated.RData")

use_condaenv('env_py3.8')
loompy <- reticulate::import('loompy')

# seurat to anndata
sceasy::convertFormat(seurat.merged, from="seurat", to="anndata",
                      outFile='SimpleMerge_Mouse_noForce_Log2_GRCh38_052424_Seurat_simpleMerge_LogNormalize_Annotated.h5ad')

# anndata to seurat
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                      outFile='filename.rds')

