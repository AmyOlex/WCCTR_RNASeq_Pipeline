##Code for Signature scores and highest Calls
#Aatish Thennavan Perou Lab
## EDITED by Amy Olex
## 6/28/23
## For subtyping Harrell PDX data.

library(Seurat)
library(optparse)

options(future.globals.maxSize = 100000 * 1024^2)

start_time <- Sys.time()


option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-c", "--pam50sigfile"), type="character", 
              help="Required. Input the PAM50 subtype signature file as a CSV.", metavar="character"),
  make_option(c("-i", "--infile"), type="character", 
              help="Path and file name of the input RData file with the merged sc samples.", metavar="character"),
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
print(paste("RData File: ", infile))
print(paste("Output Directory:" , outDir))

print("Starting Subtype Analysis...")
setwd(outDir)

# read in scsubtype gene signatures
sigdat <- read.csv(sigFile)
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

#Read in the single cell RDS object as 'Mydata'
load(inFile) ## must have a suriate object named seurat.merged saved
Mydata <- ScaleData(seurat.merged, features=temp_allgenes)
tocalc<-as.data.frame(Mydata@assays$RNA@scale.data)

#calculate mean scsubtype scores
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

final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

##Obtaining the highest call
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames
##get barcode column for Loupe Annotation output
finalm.sweep.t$barcode <- row.names(finalm.sweep.t)

##Writing out output files (rownames remain the same for both)
write.table(finalm.sweep.t[,c("barcode","SCSubtypeCall")], paste0(outDir,runID,"_SCSubtypes_Annotation.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(finalm.sweep.t[,c("barcode","SCSubtypeCall","Basal_SC","Her2E_SC","LumA_SC","LumB_SC")], paste0(outDir,runID,"_SCSubtypes_Scores.txt"), sep="\t", quote = FALSE, row.names = FALSE)

print("Subtype analysis completed!")
