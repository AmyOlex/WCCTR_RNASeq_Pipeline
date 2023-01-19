## Amy Olex
## 12.16.22
## Script to run the code for scTyper on the visvader datasets.

library(Seurat)
library(scTyper)
library(optparse)

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-f", "--seuratfile"), type="character", 
              help="Seurat object saved as an RData file in a variable named seurat.merged.", metavar="character"),
  
  make_option(c("-s", "--sigfile"), type="character", 
              help="A TSV file with the gene signatures to use for cell type classification.", metavar="character"),
  make_option(c("-p", "--phenofile"), type="character", 
              help="A CSV file with the sample annotations must be provided.", metavar="character"),
  
  make_option(c("-t", "--typelevel"), type="character", 
              help="Optional. The level of cell typing. Options are cell (default) or cluster.", 
              default = "cell", metavar="character"),
  
  
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

if (is.null(opt$sigfile)){
  print_help(opt_parser)
  stop("A TSV formatted gene siguature file must be provided.", call.=FALSE)
}

if (is.null(opt$phenofile)){
  print_help(opt_parser)
  stop("A CSV formatted sample phenotype file must be provided.", call.=FALSE)
}


runID <- opt$runid  ## "100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize.scTyper.Average.output.121922"
wd <- opt$outdir  ## "/vcu_gpfs2/home/harrell_lab/collaboration_scRNAseq/visvader_results_2021/scTyper_analysis/"
seuratfile <- opt$seuratfile ## "100422_Visvader-plussss-PDX_MASTER_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
sigfile <- opt$sigfile  ## "scTyperSignatures_12.19.22.tsv"
phenofile <- opt$phenofile  ## "SeuratSimpleMerge_Visvader-plus-PDX_noUNC_GRCh38_100422_scTyper.MetaData.csv"
typelevel <- opt$typelevel  ## "cell"

setwd(wd)

#load in data set
load(seuratfile)

#add a sample.name annotation
seurat.merged$sample.name <- seurat.merged$orig.ident

#load in the list of signatures to use
sigs <- read.delim(sigfile, header=TRUE)

#run scTyper
celltyped.seurat=scTyper(seurat.object=seurat.merged,
                         marker=sigs$Identifier,
                         wd = wd,
                         output.name = runID,
                         pheno.fn = phenofile,
                         qc = FALSE,
                         run.cellranger=FALSE,
                         norm.seurat=FALSE,
                         cell.typing.method="Average",
                         level=typelevel,
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
write.csv(data.frame(barcode = names(celltyped.seurat@active.ident), scTyperAvg.cellType = celltyped.seurat@meta.data$cell.type), file=paste0(runID,".csv", quote = FALSE, row.names = FALSE))

#save seurat object
saveRDS(celltyped.seurat, file = paste0(runID, "_scTyperTyped.rds"))


