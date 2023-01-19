## Amy Olex
## 12.16.22
## Script to run the code for scTyper on the visvader datasets.

library(Seurat)
library(scTyper)

setwd("/vcu_gpfs2/home/harrell_lab/collaboration_scRNAseq/visvader_results_2021/scTyper_analysis/")

#load in data set
load("100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData")

#add a sample.name annotation
seurat.merged$sample.name <- seurat.merged$orig.ident

#load in the list of signatures to use
sigs <- read.delim("scTyperSignatures_12.19.22.tsv", header=TRUE)

#run scTyper
celltyped.seurat=scTyper(seurat.object=seurat.merged,
                         marker=sigs$Identifier,
                         wd = getwd(),
                         output.name = "100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize.scTyper.Average.output.121922",
                         pheno.fn = "SeuratSimpleMerge_Visvader-plus-PDX_noUNC_GRCh38_100422_scTyper.MetaData.csv",
                         qc = FALSE,
                         run.cellranger=FALSE,
                         norm.seurat=FALSE,
                         cell.typing.method="Average",
                         level="cell",
                         run.inferCNV=FALSE,
                         proj.name = "scTyper",
                         gene.ref.gtf="/vcu_gpfs2/home/mccbnfolab/ref_genomes/CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf",
                         feature.to.test = "tissue.type",
                         cells.test_excluded=c("Normal"),
                         cells.test_reference = "Normal",
                         malignant.cell.type="TN",
                         report.mode=FALSE,
                         mc.cores = 30)

#save cell.type annotation to CSV file
write.csv(data.frame(barcode = names(celltyped.seurat@active.ident), scTyperAvg.cellType = celltyped.seurat@meta.data$cell.type), file="100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_scTyperAverage.cellType_12.19.22.csv", quote = FALSE, row.names = FALSE)

#save seurat object
saveRDS(celltyped.seurat, file = "100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_scTyperTyped_12.19.22.rds")


