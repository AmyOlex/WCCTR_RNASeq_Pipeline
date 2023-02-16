## Amy Olex
## 12.20.22
## inferCNV cluster script
## https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
## 

library(infercnv)
library(limma)
library(Biobase)
library(optparse)
library(Seurat)

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  
  make_option(c("-i", "--seuratfile"), type="character", default=NULL,
              help="Seurat object saved as an RData file in a variable named seurat.merged.", metavar="character"),
  
  make_option(c("-k", "--subsetK"), type="character", default=-1,
              help="Number of samples to include per group. Default is -1, which uses all samples.", metavar="integer"),
  
  make_option(c("-f", "--phenofile"), type="character", default=NULL,
              help="A CSV file with the sample annotations must be provided.", metavar="character"),
  
  make_option(c("-g", "--generef"), type="character", default="./genes.gtf",
              help="Location of the genes.gtf reference genome file.", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$seuratfile)){
  print_help(opt_parser)
  stop("A Seurat object saved as a RData file must be provided.", call.=FALSE)
}


if (is.null(opt$phenofile)){
  print_help(opt_parser)
  stop("A CSV formatted sample phenotype file must be provided.", call.=FALSE)
}


runID <- opt$runid  ## "inferCNV_TEST_Visvader-MiniSubset_012023"
wd <- opt$outdir  ## "/Users/alolex/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/SingleCell/debug_files"
seuratfile <- opt$seuratfile ## "100422_Visvader-MiniSubset_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
phenofile <- opt$phenofile  ## "inferCNV_TEST_Visvader-MiniSubset_012023_metaData.csv"
generef <- opt$generef
subsetK <- opt$subsetK ## 3


#runID <-  "inferCNV_TEST_Visvader-MiniSubset_012023"
#wd <-  "/Users/alolex/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/SingleCell/debug_files"
#seuratfile <-  "100422_Visvader-MiniSubset_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
#phenofile <-  "inferCNV_TEST_Visvader-MiniSubset_012023_metaData.csv"
#generef <- "genes.gtf"
#subsetK <- 3

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("Output Dir:", wd))
print(paste("Seurat File: ", seuratfile))
print(paste("Phenotype File: ", phenofile))
print(paste("Gene GTF: ", generef))
print(paste("Subset K: ", subsetK))
setwd(wd)

## Load in data
load(seuratfile)
samples <- read.csv(phenofile)
assay = "RNA"
print("data loaded...")
print(paste("Using assay:", assay))

## create new annotation for each cell
seurat.merged@meta.data$cell.group <- samples[match(seurat.merged@meta.data$orig.ident, samples[,"Sample_ID"]),"inferCNV.normal"]

############## Need to automate: 
## Subset to chosen samples and SELECTED_Normal samples
reference_list <- unique(seurat.merged$orig.ident[seurat.merged$cell.group %in% c("REFERENCE")])
sample_list <- unique(seurat.merged$orig.ident[!(seurat.merged$cell.group %in% c("REFERENCE"))])
sample_groups <- ifelse(subsetK < 0, list(sample_list), split(sample_list, ceiling(seq_along(sample_list)/subsetK)))
Idents(seurat.merged) <- "orig.ident"

## create an empty DF here to hold the pData information.
pData_combined <- data.frame()
refData <- ""

for(sg in sample_groups){
  print("Processing sample group: " )
  sample_subset <- c(sg, reference_list)
  print(sample_subset)
  
  ## subset Seurat object
  seurat.subset <- subset(seurat.merged, idents = sample_subset)
  
  ## Format gene expression input for inferCNV
  ## 1) create raw counts matrix
  raw.mat=as.matrix(seurat.subset[[assay]]@counts)
  
  ## 2) format annotation data
  annot=data.frame(cell.names=colnames(seurat.subset), group.st=seurat.subset@meta.data$cell.group)
  
  ## 3) get gene feature data (code from scTyper)
  ##making fdata object
  gene.ref.gtf=generef
  gtf=read.delim(gene.ref.gtf, comment.char = "#", header = F)
  gene.gtf=gtf[gtf$V3=="gene",]
  gene.annot=strsplit2(gene.gtf$V9, split = ";| ")
  gene.annot=data.frame("gene_id"=gene.annot[,2], "gene_name"=gene.annot[,8], "gene_source"=gene.annot[,11], "gene_biotype"=gene.annot[,14])
  head(gene.annot)
  mat.st=match(rownames(seurat.subset), gene.annot$gene_name)
  mat.st[is.na(mat.st)]=match(sub("\\.1", "", rownames(seurat.subset)[is.na(mat.st)]), gene.annot$gene_name)
  
  fdata=data.frame(gene.gtf[mat.st,c(1,4,5,7)], gene.annot[mat.st,])
  colnames(fdata)[1:4]=c("chr", "str", "end", "strand")
  rownames(fdata)=rownames(seurat.subset)
  
  ## adding gene annotations to Seurat obj
  seurat.subset[[assay]]@meta.features <- data.frame(seurat.subset[[assay]]@meta.features, fdata)
  names(seurat.subset[[assay]]@meta.features)
  gene.od=data.frame(gene=rownames(seurat.subset[[assay]]@meta.features), seurat.subset[[assay]]@meta.features[,c("chr", "str", "end")])
  fil.st=is.element(gene.od$chr, c(1:22, "X"))
  
  ## filtering to only primary chromosomes.
  raw.mat=raw.mat[fil.st,]
  gene.od=gene.od[fil.st,]
  gene.od$chr=factor(as.character(gene.od$chr), levels=c(1:22,"X"))
  od.st=order(gene.od$chr, gene.od$str)
  match(c(1:22, "X"), gene.od[od.st,]$chr)
  
  ## ###########
  ## Run inferCNV!
  ###########
  print("creating inferCNV object...")
  cells.test_reference = c("SELECTED_Normal")
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
  
  save(fil.cset, file=file.path("./", 'fil.cset.rda'))
  
  ######
  ## Create Z-score version of data minus the normal samples
  ##z_scores <- (data-mean(data))/sd(data)
  print("Creating ZScores..." )
  pData(fil.cset)$cnv.Zscore <- (pData(fil.cset)$cnv.score - mean(pData(fil.cset)$cnv.score))/sd(pData(fil.cset)$cnv.score)
  pData(fil.cset)$cnv.ZscoreCat <- round(pData(fil.cset)$cnv.Zscore, 0)
  pData_noref <- pData(fil.cset)[pData(fil.cset)$orig.ident %in% sg,]
  
  ##merge sample DFs together
  pData_combined <- rbind(pData_combined, pData_noref)
  
  ##save most recent ref data values (maybe average these one day?)
  refData <- pData(fil.cset)[!(pData(fil.cset)$orig.ident %in% sg),]

  
}

### Combine the pData and refData and then re-order to match the original order so we can add it as an annotation in the Seruat obj
print("Saving final files..." )

multiData <- rbind(refData, pData_combined)
multiData <- multiData[row.names(pData(fil.cset)),]

### Add to Seurat obj
seurat.merged$cnv.score=multiData$cnv.score
seurat.merged$cnv.Zscore=multiData$cnv.Zscore
seurat.merged$cnv.ZscoreCat=multiData$cnv.ZscoreCat

##write out the annotation data and files
write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvscore = seurat.merged@meta.data$cnv.score), file = paste0(celltyped.file, "_inferCNVscores.csv"), quote=FALSE, row.names=FALSE)
write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvZscore = seurat.merged@meta.data$cnv.Zscore), file = paste0(celltyped.file, "_inferCNVZscores.csv"), quote=FALSE, row.names=FALSE)
write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvZscoreCat = seurat.merged@meta.data$cnv.ZscoreCat), file = paste0(celltyped.file, "_inferCNVZscoreCat.csv"), quote=FALSE, row.names=FALSE)




# write out infercnv formatted files
#output.dir = wd
#raw.mat.fn=file.path(output.dir, paste0(runID, "_raw_counts_matrix.txt"))
#annot.fn=file.path(output.dir, paste0(runID, "_annotations.txt"))
#gene.od.fn=file.path(output.dir, paste0(runID, "_gene_order.txt"))
#write.table(raw.mat[od.st,], raw.mat.fn, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
#write.table(annot, annot.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#write.table(gene.od[od.st,], gene.od.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


###save Seurat obj with cnv scores
print("saving Seurat...")
saveRDS(seurat.merged, paste0(celltyped.file, ".inferCNVscores.rda"))

#write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvscore = seurat.merged@meta.data$cnv.score), file = paste0(celltyped.file, "_inferCNVscores.csv"), quote=FALSE, row.names=FALSE)
print("Saving images...")
png(filename = paste0(paste0(celltyped.file, ".inferCNVscores.png")))
  FeaturePlot(seurat.merged, features = "cnv.score")
dev.off()

png(filename = paste0(paste0(celltyped.file, ".inferCNVZscores.png")))
FeaturePlot(seurat.merged, features = "cnv.Zscore")
dev.off()

png(filename = paste0(paste0(celltyped.file, ".inferCNVZscoreCat.png")))
FeaturePlot(seurat.merged, features = "cnv.ZscoreCat")
dev.off()

print("COMPLETED!")
