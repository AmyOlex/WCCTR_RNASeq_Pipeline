## Amy Olex
## 12.20.22
## inferCNV cluster script
## https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
## 

library(infercnv)
library(limma)
library(Biobase)

#celltyped.file <- "100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_scTyperTyped_12.19.22.rds"
celltyped.file <- "100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
metadata.file <- "SeuratSimpleMerge_Visvader-plus-PDX_noUNC_GRCh38_100422_scTyper.MetaData.csv"
assay = "RNA"

## Load in data
load(celltyped.file)
samples <- read.csv(metadata.file)
print("data loaded...")

## create new annotation for each cell
seurat.merged@meta.data$cell.group <- samples[match(seurat.merged@meta.data$orig.ident, samples[,"Sample_ID"]),"TissueType"]

## 1) create raw counts matrix
raw.mat=as.matrix(seurat.merged[[assay]]@counts)

## 2) format annotation data
annot=data.frame(cell.names=colnames(seurat.merged), group.st=seurat.merged@meta.data$cell.group)

## 3) get gene feature data (code from scTyper)
##making fdata object
gene.ref.gtf="/vcu_gpfs2/home/mccbnfolab/ref_genomes/CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
gtf=read.delim(gene.ref.gtf, comment.char = "#", header = F)
gene.gtf=gtf[gtf$V3=="gene",]
gene.annot=strsplit2(gene.gtf$V9, split = ";| ")
gene.annot=data.frame("gene_id"=gene.annot[,2], "gene_name"=gene.annot[,8], "gene_source"=gene.annot[,11], "gene_biotype"=gene.annot[,14])
head(gene.annot)
mat.st=match(rownames(seurat.merged), gene.annot$gene_name)
mat.st[is.na(mat.st)]=match(sub("\\.1", "", rownames(seurat.merged)[is.na(mat.st)]), gene.annot$gene_name)

fdata=data.frame(gene.gtf[mat.st,c(1,4,5,7)], gene.annot[mat.st,])
colnames(fdata)[1:4]=c("chr", "str", "end", "strand")
rownames(fdata)=rownames(seurat.merged)

## adding gene annotations to Seurat obj
seurat.merged[[assay]]@meta.features <- data.frame(seurat.merged[[assay]]@meta.features, fdata)
names(seurat.merged[[assay]]@meta.features)
gene.od=data.frame(gene=rownames(seurat.merged[[assay]]@meta.features), seurat.merged[[assay]]@meta.features[,c("chr", "str", "end")])
fil.st=is.element(gene.od$chr, c(1:22, "X"))

## filtering to only primary chromosomes.
raw.mat=raw.mat[fil.st,]
gene.od=gene.od[fil.st,]
gene.od$chr=factor(as.character(gene.od$chr), levels=c(1:22,"X"))
od.st=order(gene.od$chr, gene.od$str)
match(c(1:22, "X"), gene.od[od.st,]$chr)

# write out infercnv formatted files
output.dir = "./"
raw.mat.fn=file.path(output.dir, "raw_counts_matrix.txt")
annot.fn=file.path(output.dir, "annotations.txt")
gene.od.fn=file.path(output.dir, "gene_order.txt")
write.table(raw.mat[od.st,], raw.mat.fn, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(annot, annot.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(gene.od[od.st,], gene.od.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

###########
## Run inferCNV!
###########
print("creating inferCNV object...")
cells.test_reference = c("Normal")
infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=raw.mat[od.st,],
                                              annotations_file=annot.fn,
                                              gene_order_file=gene.od.fn,
                                              ref_group_names=cells.test_reference)
## This took 3.5 hours for about 10k cells!
print("running inferCNV command...")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_infercnv",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             num_threads=30)

save(infercnv_obj, file=file.path("./", 'infercnv_obj.o.rda'))

print("completed inferCNV...")

### Filter out by standard deviations in reference groups
ref.sd=apply(infercnv_obj@expr.data[,infercnv_obj@reference_grouped_cell_indices[[1]]], 1, sd)
sd.cut = 0.3
sd.fil.st <- ref.sd < sd.cut
infercnv_obj@expr.data=infercnv_obj@expr.data[sd.fil.st,]
infercnv_obj@count.data=infercnv_obj@count.data[sd.fil.st,]
infercnv_obj@gene_order=infercnv_obj@gene_order[sd.fil.st,]

save(infercnv_obj, file=file.path("./", 'infercnv_obj.sd_filtered.rda'))
save(infercnv_obj, file=file.path("./", 'infercnv_obj.rda'))

#make filtered cnvset

fil.cset=make.eset(expr = infercnv_obj@expr.data, pdata = seurat.merged@meta.data[colnames(infercnv_obj@expr.data),], fdata = infercnv_obj@gene_order)

pData(fil.cset)$cnv.score=apply(exprs(fil.cset)-1, 2, function(a) sum(abs(a)))

seurat.merged$cnv.score=pData(fil.cset)$cnv.score
save(fil.cset, file=file.path("./", 'fil.cset.rda'))

###save Seurat obj with cnv scores
print("saving Seurat...")
saveRDS(seurat.merged, paste0(celltyped.file, ".inferCNVscores.rda"))

write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvscore = seurat.merged@meta.data$cnv.score), file = paste0(celltyped.file, "_inferCNVscores.csv"), quote=FALSE, row.names=FALSE)

png(filename = paste0(paste0(celltyped.file, ".inferCNVscores.png")))
  FeaturePlot(seurat.merged, features = "cnv.score")
dev.off()

print("COMPLETED!")
