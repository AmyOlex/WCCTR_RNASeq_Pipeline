## Code for Quantile Signature Scoring
## Code based on SCSubtype HighestCalls script by Aatish Thennavan Perou Lab
## Amy Olex
## 03/12/24
## For scoring cells based on expressed gene signatures.
## Gene signatures MUST be made up of ONLY Up-regulated/Highly Expressed Genes


library(Seurat)
library(optparse)
library(ComplexHeatmap)
library(dplyr)

options(future.globals.maxSize = 100000 * 1024^2)

start_time <- Sys.time()

debug=FALSE
if(debug){
  runID <- "032124_PyMT_Primary_Pooled_withLuc_AmbientAdj"
  #inFile <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/SignatureScoring_Testing/PyMT_Primary_Master_Merge/092623_PrimaryTumor_PyMT_Master_SimpleMerge_mm10firefly_noLuc_Seurat_simpleMerge_LogNormalize_Annotated.rds"
  inFile <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/SignatureScoring_Testing/PyMT_Primary_Master_merge_withLuc/031924_PrimaryTumor_pooled_PyMT_Master_SimpleMerge_mm10firefly_ambiantadj_withLuc_Seurat_simpleMerge_LogNormalize_Annotated.rds"
  sigFile <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/SignatureScoring_Testing/Bos_Gene_Signatures.csv"
  outDir <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/SignatureScoring_Testing/PyMT_Primary_Master_Merge_withLuc/"
  setwd(outDir)
  }



option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-c", "--sigfile"), type="character", 
              help="Required. Input the gene signature file as a CSV. One list of Gene Symbols per column with the column name as the Signature Name.", metavar="character"),
  make_option(c("-i", "--infile"), type="character", 
              help="Path and file name of the input RDS file with the Seurat object for analysis.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$pam50sigfile)){
  print_help(opt_parser)
  stop("A PAM50 signature file must be provided.", call.=FALSE)
}

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("An input RData file with Seurat single cell object must be provided.", call.=FALSE)
}


runID <- opt$runid
sigFile <- opt$pam50sigfile
outDir <- paste0(opt$outdir,"/")
inFile <- opt$infile

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("SignatureFile: ", sigFile))
print(paste("RData File: ", inFile))
print(paste("Output Directory:" , outDir))

print("Starting Subtype Analysis...")
setwd(outDir)

# read in scsubtype gene signatures
sigdat <- read.csv(sigFile)
temp_allgenes <- c(as.vector(sigdat[,"IFN.gamma"]),
                   as.vector(sigdat[,"ECM.EO771"]),
                   as.vector(sigdat[,"ECM.PyMT"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

#Read in the single cell RDS object as 'Mydata'
print("Loading RData file...")
seurat.obj <- readRDS(inFile) ## must have a Seurat object in it.

## Recalculate the nCount_RNA
seurat.obj$nCount_RNA <- colSums(seurat.obj@assays$RNA$counts) 

print("RDS file loaded, scaling data...")
Mydata <- ScaleData(seurat.obj, features=temp_allgenes, model.use = "negbinom", 
                    vars.to.regress = "nCount_RNA", center=FALSE, scale=TRUE)
tocalc<-as.data.frame(Mydata@assays$RNA@scale.data)

#calculate mean signature scores
print("Calculating mean signature scores...")
outdat <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))

outdat_sum <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))

## I tried removing any signature gene that had a sum of zero across all cells, but this is the scaled data, so it is never all zero.
## Maybe if I round the absolute values? Or say anything below0.5, but that is really arbitrary.
for(i in 1:ncol(sigdat)){
  
  siggenes <- as.character(sigdat[,i])
  siggenes<-unique(siggenes[siggenes != ""])
  genes<-which(rownames(tocalc) %in% siggenes)
  
  ## remove any genes in the signature that have no expression across all cells
  siggene_data <- tocalc[genes,]
  siggene_data <- siggene_data[rowSums(abs(siggene_data))>0,]
  
  ## Using average
  temp<-apply(siggene_data,2,function(x){mean(as.numeric(x),na.rm=FALSE)})
  outdat[i,]<-as.numeric(temp)
  
  ## Using SUM
  temp_sum<-apply(siggene_data,2,function(x){sum(as.numeric(x),na.rm=FALSE)})
  outdat_sum[i,]<-as.numeric(temp_sum)
  
}

outdat_t <- as.data.frame(t(outdat))
outdat_sum_t <- as.data.frame(t(outdat_sum))

new_names <- c(paste0(names(outdat_t), "_SigScoreAVG"), paste0(names(outdat_t), "_SigCatAVG"))
numcol <- ncol(outdat_t)

new_names_sum <- c(paste0(names(outdat_sum_t), "_SigScoreSUM"), paste0(names(outdat_sum_t), "_SigCatSUM"))
#numcol <- ncol(outdat_sum_t)



## CELL SIGNATURE CLASSIFICATION for AVG
for(i in 1:ncol(outdat_t)){
  this_quant <- quantile(outdat_t[,i])
  outdat_t[,i+numcol] <- 'LessThanQ3'
  outdat_t[outdat_t[,i]>this_quant[4],i+numcol] <- 'HigherThanQ3'
}

names(outdat_t) <- new_names
outdat_t$barcode <- row.names(outdat_t)


## CELL SIGNATURE CLASSIFICATION for SUM
for(i in 1:ncol(outdat_sum_t)){
  this_quant <- quantile(outdat_sum_t[,i])
  outdat_sum_t[,i+numcol] <- 'LessThanQ3'
  outdat_sum_t[outdat_sum_t[,i]>this_quant[4],i+numcol] <- 'HigherThanQ3'
}

names(outdat_sum_t) <- new_names_sum
outdat_sum_t$barcode <- row.names(outdat_sum_t)

## write out annotation files for Loupe and add each to the Seurat Object
for(i in 1:(ncol(outdat_t)-1)){
  write.table(file=paste0(runID, "_SigAnnotation_scaled_",new_names[i],".csv"), outdat_t[,c("barcode",new_names[i])], row.names = FALSE, quote = FALSE, sep=",")
  write.table(file=paste0(runID, "_SigAnnotation_scaled_",new_names_sum[i],".csv"), outdat_sum_t[,c("barcode",new_names_sum[i])], row.names = FALSE, quote = FALSE, sep=",")
  
  seurat.obj[[new_names[i]]] <- outdat_t[,new_names[i]]
  seurat.obj[[new_names_sum[i]]] <- outdat_sum_t[,new_names_sum[i]]
}

## Save update Seurat Object
saveRDS(seurat.obj, file = paste0(runID, "_SeuratObj_wSignatureAnnots.rds"))

## GENERATE HEATMAPS
for(i in 1:ncol(sigdat)){
  
  siggenes <- as.character(sigdat[,i])
  siggenes<-unique(siggenes[siggenes != ""])
  genes<-which(rownames(tocalc) %in% siggenes)
  
  this_sig_exp <- as.data.frame(t(tocalc[genes,]))
  this_sig_exp$classes <- outdat_t[,new_names[i+dim(outdat)[1]]]
  this_sig_exp$Row.names <- row.names(this_sig_exp)
  
  ## Subset to 20% in each group
  pop <- this_sig_exp %>% group_by(classes)
  subset <- slice_sample(pop, prop=0.2)
  row.names(subset) <- subset$Row.names
  subset <- subset[rowSums(abs(subset[,1:(ncol(subset)-2)]))>0,]
  #subset <- subset[,colSums(abs(subset[,1:(ncol(subset)-2)]))>0]
  
  subset_scaled <- na.omit(as.matrix(t(scale(subset[,1:(ncol(subset)-2)], center = TRUE, scale = FALSE))))
  #val <- max(abs(subset[,1:(ncol(subset)-2)]))
  val <- max(abs(subset_scaled))
  breaks = c((-1*val),0,val)
  
  # Create a color mapping function centered on zero
  color_mapping = circlize::colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))
  
  
  #colAnnot <- HeatmapAnnotation(Class = outdat_t[,new_names[i+dim(outdat)[1]]], col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
  colAnnot <- HeatmapAnnotation(Class = subset$classes, col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
  
  
  png(filename=paste0(runID,"_",names(sigdat)[i],"_heatmap_centered.png")) 
    print(Heatmap(subset_scaled, col=color_mapping, top_annotation = colAnnot, column_split = subset$classes)) #, heatmap_legend_param = list(Class = "Class", )))  
  dev.off()
  
  val <- max(abs(subset[,1:(ncol(subset)-2)]))
  breaks = c((-1*val),0,val)
  
  # Create a color mapping function centered on zero
  color_mapping = circlize::colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))
  
  
  png(filename=paste0(runID,"_",names(sigdat)[i],"_heatmap.png"))
    print(Heatmap(t(subset[,1:(ncol(subset)-2)]), col=color_mapping, top_annotation = colAnnot, column_split = subset$classes)) #, heatmap_legend_param = list(Class = "Class", )))  
  dev.off()
}






######################################## NO, Don't Use This Part ############################
#### Using the SUM of RAW data

## Recalculate the nCount_RNA
#Mydata$nCount_RNA <- colSums(Mydata@assays$RNA$counts)

## Remove all genes with a zero count across all cells
#Mydata_filt <- Mydata[rowSums(Mydata@assays$RNA@counts) > 0, ]

## Calculate sums
#tocalc_sum <- as.data.frame(Mydata_filt@assays$RNA@counts)

# print("Calculating SUM signature scores...")
# outdat_sum <- matrix(0,
#                  nrow=ncol(sigdat),
#                  ncol=ncol(tocalc_sum),
#                  dimnames=list(colnames(sigdat),
#                                colnames(tocalc_sum)))
# 
# outdat_sum_norm <- outdat_sum
# 
# for(i in 1:ncol(sigdat)){
# 
#   siggenes <- as.character(sigdat[,i])
#   siggenes <- unique(siggenes[siggenes != ""])
#   genes<-which(rownames(tocalc_sum) %in% siggenes)
# 
#   temp<-apply(tocalc_sum[genes,],2,function(x){sum(as.numeric(x),na.rm=TRUE)})
# 
#   ## I am thinking I should probably normalize by the nCount_RNA, which is the UMI basically for read depth of each cells.
#   outdat_sum[i,]<-as.numeric(temp)
#   outdat_sum_norm[i,]<-as.numeric(temp)/Mydata_filt$nCount_RNA
# 
# }
# 
# outdat_sum_t <- as.data.frame(t(outdat_sum))
# outdat_sum_norm_t <- as.data.frame(t(outdat_sum_norm))
# 
# ### Assess gene expression for EO771 signature
# #siggenes <- as.character(sigdat[,2])
# #siggenes<-unique(siggenes[siggenes != ""])
# #genes<-which(rownames(tocalc_sum) %in% siggenes)
# 
# #temp<-rowSums(tocalc_sum[genes,])
# #View(as.data.frame(temp))
# 
# ## CELL SIGNATURE CLASSIFICATION
# for(i in 1:ncol(outdat_sum_norm_t)){
#   this_quant <- quantile(outdat_sum_norm_t[,i])
#   outdat_sum_norm_t[,i+numcol] <- 'LessThanQ3'
#   outdat_sum_norm_t[outdat_sum_norm_t[,i]>this_quant[4],i+numcol] <- 'HigherThanQ3'
# }
# 
# names(outdat_sum_norm_t) <- new_names
# outdat_sum_norm_t$barcode <- row.names(outdat_sum_norm_t)
# 
# ## write out annotation files for Loupe
# for(i in (1+dim(outdat_sum)[1]):(ncol(outdat_sum_norm_t)-1)){
#   write.table(file=paste0(runID, "_NormSum_SigCategories_scaled_",new_names[i],".csv"), outdat_sum_norm_t[,c("barcode",new_names[i])], row.names = FALSE, quote = FALSE, sep=",")
# }
# 
# 
# ## Generate heatmaps
#   ###  This heat map is supposed to show the raw count data from the single cell, but I keep getting Nan's
#   ###  Maybe need to re-evaluate the Seurat Scale method and see if I can get it to center the data.
# for(i in 1:ncol(sigdat)){
# 
#   siggenes <- as.character(sigdat[,i])
#   siggenes<-unique(siggenes[siggenes != ""])
#   genes<-which(rownames(tocalc_sum) %in% siggenes)
# 
#   tmp <- tocalc_sum[genes,]
#   tmp_filt <- tmp[rowSums(tmp)!=0,]
# 
#   this_sig_exp <- as.data.frame(t(tmp_filt))
#   this_sig_exp$classes <- outdat_sum_norm_t[,new_names[i+dim(outdat_sum)[1]]]
#   this_sig_exp$Row.names <- row.names(this_sig_exp)
# 
#   ## Subset to 20% in each group
#   pop <- this_sig_exp %>% group_by(classes)
#   subset <- slice_sample(pop, prop=0.2)
#   row.names(subset) <- subset$Row.names
#   subset <- subset[rowSums(subset[,1:(ncol(subset)-2)])>0,]
# 
#   subset_scaled <- na.omit(as.matrix(t(scale(subset[,1:(ncol(subset)-2)], center = TRUE, scale = TRUE))))
#   val <- max(abs(subset_scaled))
#   breaks = c((-1*val),0,val)
# 
#   # Create a color mapping function centered on zero
#   color_mapping = circlize::colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))
# 
# 
#   #subset_scaled[is.na(subset_scaled)] <- 0
#   #colAnnot <- HeatmapAnnotation(Class = outdat_t[,new_names[i+dim(outdat)[1]]], col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
#   colAnnot <- HeatmapAnnotation(Class = subset$classes, col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
# 
#   #print(Heatmap(scale(t(subset[,1:(ncol(subset)-2)]), center = TRUE, scale = FALSE), top_annotation = colAnnot, column_split = subset$classes))
# 
#   png(filename=paste0(runID,"_",names(sigdat)[i],"_NormSum_heatmap_centered.png"))
#     print(Heatmap(subset_scaled, col = color_mapping, top_annotation = colAnnot, column_split = subset$classes)) #, heatmap_legend_param = list(Class = "Class", )))
#   dev.off()
# 
#   png(filename=paste0(runID,"_",names(sigdat)[i],"_NormSum_heatmap.png"))
#     print(Heatmap(t(subset[,1:(ncol(subset)-2)]),  col=c("white", "red"), top_annotation = colAnnot, column_split = subset$classes)) #, heatmap_legend_param = list(Class = "Class", )))
#   dev.off()
# }

# 
# 
# 
# ########################### Not using the code below yet ####################
# 
# 
# #final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
# #final<-as.data.frame(final)
# #is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
# #finalm<-as.matrix(final)
# 
# ##Scaling scores function before calling the highest Call
# print("Scaling scores function before calling the highest Call...")
# center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
#   get_average <- function(v) sum(v * row.w)/sum(row.w)
#   average <- apply(x, 2, get_average)
#   sweep(x, 2, average)
# }
# 
# print("Obtaining the highest call...")
# finalmt<-as.data.frame(t(finalm))
# finalm.sweep.t<-center_sweep(finalmt)
# Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
# finalm.sweep.t$SCSubtypeCall <- Finalnames
# ##get barcode column for Loupe Annotation output
# finalm.sweep.t$barcode <- row.names(finalm.sweep.t)
# 
# print("Saving data...")
# ##Writing out output files (rownames remain the same for both)
# write.table(finalm.sweep.t[,c("barcode","SCSubtypeCall")], paste0(outDir,runID,"_SCSubtypes_Annotation.txt"), sep="\t", quote = FALSE, row.names = FALSE)
# write.table(finalm.sweep.t[,c("barcode","SCSubtypeCall","Basal_SC","Her2E_SC","LumA_SC","LumB_SC")], paste0(outDir,runID,"_SCSubtypes_Scores.txt"), sep="\t", quote = FALSE, row.names = FALSE)
# 
# print("Subtype analysis completed!")
