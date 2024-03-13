## Code for Quantile Signature Scoring
## Code based on SCSubtype HighestCalls script by Aatish Thennavan Perou Lab
## Amy Olex
## 03/12/24
## For scoring cells based on expressed gene signatures.
## Gene signatures MUST be made up of ONLY Up-regulated/Highly Expressed Genes


library(Seurat)
library(optparse)

options(future.globals.maxSize = 100000 * 1024^2)

start_time <- Sys.time()

debug=FALSE
if(debug){
  runID <- "sigscoreTEST"
  inFile <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/CellChat/091923_BrainMetMaster_SimpleMerge_mm10firefly_noLuc_Seurat_simpleMerge_LogNormalize_Annotated.rds"
  sigFile <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/SignatureScoring_Testing/Bos_Gene_Signatures.csv"
  outDir <- "/Users/alolex/Desktop/CCTR_LOCAL_Analysis_noBackups/PaulaBos_LocalWorkDir/SignatureScoring_Testing/"
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

#print("RData file loaded, scaling data...")
Mydata <- ScaleData(seurat.obj, features=temp_allgenes)
tocalc<-as.data.frame(Mydata@assays$RNA@scale.data)

print("RData file loaded, skipping scaling and using raw RNA reads...")
#tocalc<-as.data.frame(seurat.obj@assays$RNA@counts)

#calculate mean scsubtype scores
print("Calculating mean signature scores...")
outdat <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))
for(i in 1:ncol(sigdat)){
  
  # sigdat[i,!is.na(sigdat[i,])]->module
  # row <- as.character(unlist(module))
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
  
  outdat[i,]<-as.numeric(temp)
  
}

outdat_t <- as.data.frame(t(outdat))
## Score scaling
## NO, I don't think we should scale.  I need to see the heatmap.
##outdat_t <-as.data.frame(t(scale(outdat, center = TRUE, scale = FALSE)))

new_names <- c(names(outdat_t), paste0(names(outdat_t), "_SigScore"))
numcol <- ncol(outdat_t)




## CELL SIGNATURE CLASSIFICATION
for(i in 1:ncol(outdat_t)){
  this_quant <- quantile(outdat_t[,i])
  outdat_t[,i+numcol] <- 'LessThanQ3'
  outdat_t[outdat_t[,i]>this_quant[4],i+numcol] <- 'HigherThanQ3'
}

names(outdat_t) <- new_names
outdat_t$barcode <- row.names(outdat_t)

## write out annotation files for Loupe
for(i in (1+dim(outdat)[1]):(ncol(outdat_t)-1)){
  write.table(file=paste0("SigCategories_scaled_",new_names[i],".csv"), outdat_t[,c("barcode",new_names[i])], row.names = FALSE, quote = FALSE, sep=",")
}

library(ComplexHeatmap)
library(dplyr)
## Generate heatmaps
for(i in 1:ncol(sigdat)){
  
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  
  this_sig_exp <- as.data.frame(t(tocalc[genes,]))
  this_sig_exp$classes <- outdat_t[,new_names[i+dim(outdat)[1]]]
  this_sig_exp$Row.names <- row.names(this_sig_exp)
  
  ## Subset to 10% in each group
  pop <- this_sig_exp %>% group_by(classes)
  subset <- slice_sample(pop, prop=0.01)
  row.names(subset) <- subset$Row.names

  #colAnnot <- HeatmapAnnotation(Class = outdat_t[,new_names[i+dim(outdat)[1]]], col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
  colAnnot <- HeatmapAnnotation(Class = subset$classes, col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
  
  Heatmap(t(subset[,1:(ncol(subset)-2)]), top_annotation = colAnnot)  
  
}



#### Using the SUM of RAW data
tocalc_sum <- as.data.frame(Mydata@assays$RNA@counts)

print("Calculating SUM signature scores...")
outdat_sum <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc_sum),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc_sum)))
for(i in 1:ncol(sigdat)){
  
  # sigdat[i,!is.na(sigdat[i,])]->module
  # row <- as.character(unlist(module))
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc_sum) %in% row)
  
  temp<-apply(tocalc_sum[genes,],2,function(x){sum(as.numeric(x),na.rm=TRUE)})
  
  outdat_sum[i,]<-as.numeric(temp)
  
}

outdat_sum_t <- as.data.frame(t(outdat_sum))

### Assess gene expression for EO771 signature
row <- as.character(sigdat[,2])
row<-unique(row[row != ""])
genes<-which(rownames(tocalc_sum) %in% row)

temp<-rowSums(tocalc_sum[genes,])
View(as.data.frame(temp))


temp2 <- 

## Generate heatmaps
  ###  This heat map is supposed to show the scaled raw count data from the single cell, but I keep getting Nan's 
  ###  Maybe need to re-evaluate the Seurat Scale method and see if I can get it to center the data.
for(i in 1:ncol(sigdat)){
  
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc_sum) %in% row)
  
  tmp <- tocalc_sum[genes,]
  tmp_filt <- tmp[rowSums(tmp)!=0,]
  
  this_sig_exp <- as.data.frame(t(tmp_filt))
  this_sig_exp$classes <- outdat_t[,new_names[i+dim(outdat)[1]]]
  this_sig_exp$Row.names <- row.names(this_sig_exp)
  
  ## Subset to 10% in each group
  pop <- this_sig_exp %>% group_by(classes)
  subset <- slice_sample(pop, prop=0.01)
  row.names(subset) <- subset$Row.names
  subset_scaled <- as.matrix(scale(t(subset[,1:(ncol(subset)-2)]), center = TRUE, scale = TRUE))
  
  #colAnnot <- HeatmapAnnotation(Class = outdat_t[,new_names[i+dim(outdat)[1]]], col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
  colAnnot <- HeatmapAnnotation(Class = subset$classes, col = list(Class = c("HigherThanQ3" = "red", "LessThanQ3" = "blue")))
  
  Heatmap(subset_scaled, top_annotation = colAnnot)  
  
}




########################### Not using the code below yet ####################


#final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
#final<-as.data.frame(final)
#is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
#finalm<-as.matrix(final)

##Scaling scores function before calling the highest Call
print("Scaling scores function before calling the highest Call...")
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

print("Obtaining the highest call...")
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames
##get barcode column for Loupe Annotation output
finalm.sweep.t$barcode <- row.names(finalm.sweep.t)

print("Saving data...")
##Writing out output files (rownames remain the same for both)
write.table(finalm.sweep.t[,c("barcode","SCSubtypeCall")], paste0(outDir,runID,"_SCSubtypes_Annotation.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(finalm.sweep.t[,c("barcode","SCSubtypeCall","Basal_SC","Her2E_SC","LumA_SC","LumB_SC")], paste0(outDir,runID,"_SCSubtypes_Scores.txt"), sep="\t", quote = FALSE, row.names = FALSE)

print("Subtype analysis completed!")
