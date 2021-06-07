## Amy Olex
## 4/5/2021
## File to merge dead cell filtered scRNASeq data from Paula's experiments
## Steps for Merging:
## 1) create Seurat .h5Seurat files for each sample
## 2) Rename Barcodes
## 3) Normalize each data set using SCT or Log
## 4) Create a list of Seurat objects to merge
## 5) Select integration features and find anchors
## 6) Merge Data Sets
## 7) Write out 10X files
## 8) Create sample annotations
## 9) Run PCA, tSNE, and UMAP, and output feature annotations
## 10) Cluster cells and output feature annotations
## 11) Save Merged and annotated Seurat Object as an h5Seurat file.

## INPUTS
## Path to deadcell filtered h5Seurat file OR path to 10X filtered_features_bc_matrix folder
## Flag to delineate input types.
## 
## sample file format:
## SampleName, DataType (Seurat|10X), SamplePath, Condition

## UPDATED 6/4/21 to also provide a simple merge option

library(Seurat)
#library(SeuratDisk)
library(stringr)
library(DropletUtils)
library(optparse)
library(foreach)
library(doParallel)
library(future)

localtest = FALSE
###########################################
#### Local Testing Block
if(localtest){
  runID <- "TEST"
  inFile <- "./configtest.csv"
  outDir <- "./"
  savedir <- paste0(outDir,runID)
  numCores <- 2
  numAnchors <- 2000
  normalization <- "LogNormalize"
  mergeType <- "simple"
  parallel <- FALSE
  saveH5 <- TRUE
  options(future.globals.maxSize = 3000 * 1024^2)
}
###########################################

#numCores <- detectCores()
#registerDoParallel(numCores)

options(future.globals.maxSize = 10000 * 1024^2)

start_time <- Sys.time()


option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-c", "--configfile"), type="character", 
              help="Required. Input config file with sample information.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-i", "--normalization"), type="character", 
              help="Type of normalization to run prior to merging/integrating (SCT or LogNormalize). Default is SCT.",
              default = "SCT", metavar="character"),
  make_option(c("-t", "--type"), type="character", 
              help="Type of merging to do (simple or integration). Default is simple.",
              default = "simple", metavar="character"),
  make_option(c("--parallel"), type="logical", 
              help="Optional. Use paralellization for Merging and Integration (note, normalization of individual samples is always parallelized).",
              default = FALSE, action = "store_true", metavar="logical"),
  make_option(c("-n", "--numCores"), type="numeric", 
              help="Optional. Specify number of cores to use. Default is 1.",
              default = 1, metavar="character"),
  make_option(c("-a", "--numAnchors"), type="numeric", 
              help="Optional. Specify number of anchor genes to use. Default is 2000.",
              default = 2000, metavar="character"),
  make_option(c("--saveH5"), type="logical", 
              help="Optional. Save merged and annotated Seruat object as a .h5Seruat file. Default saves as a .RData file.",
              default = FALSE, action = "store_true", metavar="logical")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$configfile)){
  print_help(opt_parser)
  stop("A configuration file with sample information must be provided.", call.=FALSE)
}


runID <- opt$runid
inFile <- opt$configfile
outDir <- opt$outdir
savedir <- paste0(outDir,runID)
numCores <- opt$numCores
numAnchors <- opt$numAnchors
normalization <- opt$normalization
parallel <- opt$parallel
saveH5 <- opt$saveH5
mergeType <- opt$type

#### Ok, print out all the current running options as a summary.

print("Summary of input options:\n")
print(paste("Run Name:", runID))
print(paste("ConfigFile:", inFile))
print(paste("Merge Type:", mergeType))
print(paste("Output Directory:", outDir))
print(paste("Cores Specified:", numCores))
print(paste("Number of Anchor Genes:", numAnchors))
print(paste("Integrantion/Normalization Type:", normalization))
print(paste("Using Parallel: ", parallel))
print(paste("Saving h5Seurat file:", saveH5))

if(!parallel){
  print("WARNING: Seurat merge and integration functions will NOT be parallized.  Use the --parallel flag to parallelize these functions.")
}

############## PROCESS CONFIG FILE###########################
# Read in the provided config file and loop for each row.
toProcess = read.table(inFile, header=TRUE, sep=",", stringsAsFactors = FALSE)
#toProcess = read.table("/Users/alolex/Desktop/CCTR_Git_Repos/PBos_scRNASeq/MergeScriptTesting/configtest.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
toProcess$SeuratObj <- ""
print(paste(dim(toProcess)[1], " rows were found."))

#Register num cores for doParallel loop
registerDoParallel(numCores)

sc_start <- Sys.time()

seurat_list <- foreach(i=1:dim(toProcess)[1]) %dopar% {
  print(paste("Importing data for row", i, "from sample", toProcess[i,"SampleName"]))
  
  if(toProcess[i,"DataType"] == "Seurat"){
    print("Loading Seurat h5 File...")
    library(SeuratDisk)
    h5 <- LoadH5Seurat(toProcess[i,"SamplePath"])
  }
  else if(toProcess[i,"DataType"] == "10X"){
    print("Loading 10X feature matrix...")
    # Load the data set and create the Seurat object
    sc.data <- Read10X(data.dir = toProcess[i,"SamplePath"])  ## must point to the filtered_feature_bc_matrix directory.
    # Initialize the Seurat object
    h5 <- CreateSeuratObject(counts = sc.data, project = toProcess[i,"SampleName"], min.cells = 0, min.features = 0)
    
  }
  else{
    print("Error unknown data type.")
    quit(1)
  }
  
  print("Renaming Cells...")
  ## Rename the barcodes in each file with a count greater than 1
  if(i > 1){
    h5 <- RenameCells(h5,  new.names = str_replace(names(h5$orig.ident), "-1", paste0("-",i)))
  }
  
  mid_time <- Sys.time()
  
  if(normalization == "SCT"){
    print("SCTransform...")
    ## Normalize each dataset with SCT
    h5 <- SCTransform(h5, verbose = FALSE)
  }
  else if(normalization == "LogNormalize"){
    print("LogNormalize...")
    h5 <- NormalizeData(h5, normalization.method = "LogNormalize", verbose = FALSE)
    h5 <- FindVariableFeatures(h5, selection.method = "vst", nfeatures = numAnchors, verbose = FALSE)
    h5 <- ScaleData(h5, verbose = FALSE)
  }
  
  
  ## Add seurat obj to list now
  #seurat_list <- c(seurat_list,h5)  
  
  print(Sys.time() - mid_time)
  
  h5
}  ## end processing of each sample.

print("Total Time to run sample normalization:")
print(Sys.time() - sc_start)

print("Sample import and normalization completed...")
print("Finding integration features...")
## Find integration features

#seurat_list
if(parallel){
  ## Set Future Plan for asynchronous execution
  plan("multiprocess", workers = numCores)
}

plan()


###### Run Integration
if(mergeType == "integration"){
  print("Running Integration...")
  ## Select integration features
  if(normalization == "SCT"){
    mid_time <- Sys.time()
    hfile_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = numAnchors)
    
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = hfile_features, verbose = FALSE)
    print(Sys.time() - mid_time)
  } else{
    hfile_features <- numAnchors
  }
  ## Find the anchors and then integrate the data sets.
  mid_time <- Sys.time()
  print("Finding Anchors...")
  #seurat.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = integrationNormalization, anchor.features = hfile_features, reference = 5, verbose = FALSE)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = integrationNormalization, anchor.features = hfile_features, verbose = FALSE)
  
  ## Save Anchor Feature File
  write.csv(seurat.anchors@anchor.features, file=paste0(savedir, "_Seurat", normalization, "Merge_AnchorList.csv"), quote = FALSE, row.names = FALSE)
  
  print(Sys.time() - mid_time)
  
  mid_time <- Sys.time()
  print("Performing integration...")
  seurat.merged <- IntegrateData(anchorset = seurat.anchors, normalization.method = normalization, verbose = FALSE)
  print(Sys.time() - mid_time)
  
} else if(mergeType == "simple"){
  print("Running simple merge...")
  
  seurat.merged <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], merge.data=TRUE, project=runID)
  seurat.merged <- FindVariableFeatures(seurat.merged, selection.method = "vst", nfeatures = numAnchors, verbose = FALSE)
}



mid_time <- Sys.time()
print("Saving to 10X...")
## Saving to 10X format:
if(normalization == "SCT"){
  write10xCounts(x=seurat.merged@assays$SCT@data, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_SCTdata.h5"), version="3")
  write10xCounts(x=seurat.merged@assays$RNA@data, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_RNAdata.h5"), version="3")
  print(Sys.time() - mid_time)
} else if(normalization == "LogNormalize"){
  write10xCounts(x=seurat.merged@assays$RNA@data, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_RNAdata.h5"), version="3")
  print(Sys.time() - mid_time)
}

print("Adding Condition...")
mid_time <- Sys.time()
## Add Sample annotations
idents <- data.frame(barcode = names(seurat.merged@active.ident), LibraryID = seurat.merged@active.ident)
condition <- idents
names(condition) <- c("barcode","Condition")
#condition$Condition <- as.character(condition$LibraryID)

for(i in dim(toProcess)[1]){
  condition[condition$Condition == toProcess[i,1],"Condition"] <- toProcess[i,"Condition"]
}

seurat.merged$Condition <- condition[,"Condition",drop=FALSE]

## Save Sample Annotations
write.csv(idents, file=paste0(savedir, "_Seurat_",mergeType,"Merge_", normalization, "_LibraryID.csv"), quote = FALSE, row.names = FALSE)
write.csv(condition, file=paste0(savedir, "_Seurat_",mergeType,"Merge_", normalization, "_Condition.csv"), quote = FALSE, row.names = FALSE)
print(Sys.time() - mid_time)

print("Running PCA, UMAP, tSNE...")
mid_time <- Sys.time()
### Create Visualizations
if(mergeType == "integrated"){
  DefaultAssay(seurat.merged) <- "integrated"
} else {
  DefaultAssay(seurat.merged) <- "RNA"
}


if(normalization == "LogNormalize"){
  seurat.merged <- ScaleData(seurat.merged, verbose = FALSE)
}

seurat.merged <- RunPCA(seurat.merged, verbose = FALSE)
seurat.merged <- RunUMAP(seurat.merged, dims = 1:30)
seurat.merged <- RunTSNE(seurat.merged, dims = 1:30)

png(filename = "umap_test.png")
DimPlot(seurat.merged, reduction = "pca", group.by="orig.ident")
dev.off()

png(filename = "tsne_test.png")
DimPlot(seurat.merged, reduction = "tsne", group.by="orig.ident")
dev.off()

write.csv(seurat.merged@reductions$umap@cell.embeddings, file = paste0(savedir, "_UMAPCoordinates_30PCs_",mergeType,"Merge_",normalization,".csv"), quote = FALSE)
write.csv(seurat.merged@reductions$tsne@cell.embeddings, file = paste0(savedir, "_tSNECoordinates_30PCs_",mergeType,"Merge_",normalization,".csv"), quote = FALSE)
print(Sys.time() - mid_time)

print("SNN Clustering...")
mid_time <- Sys.time()
### Cluster the Data
seurat.merged <- FindNeighbors(seurat.merged, reduction = "pca", dims = 1:30, graph.name = "merged_snn")
seurat.merged <- FindClusters(seurat.merged, resolution = 0.4, graph.name = "merged_snn")

clusters <- data.frame("barcode" = names(seurat.merged$merged_snn_res.0.4), "SNN_res0.4_Clusters" = seurat.merged$merged_snn_res.0.4)
write.csv(clusters, file = paste0(savedir, "_SNN_Clusters_res0.4_",mergeType,"MergedData.csv"), quote = FALSE)
print(Sys.time() - mid_time)

print("Saving Annotated Seurat File...")
mid_time <- Sys.time()

### Save Seurat Object for future use
if(saveH5){
  library(SeuratDisk)
  SaveH5Seurat(seurat.merged, filename=paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated.h5Seurat"), overwrite = TRUE)
}else{
  save(seurat.merged, file = paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated.RData"), compress = TRUE)
}


print(Sys.time() - mid_time)

print("Seurat merging completed in: ")
end_time <- Sys.time()
print(end_time - start_time)

