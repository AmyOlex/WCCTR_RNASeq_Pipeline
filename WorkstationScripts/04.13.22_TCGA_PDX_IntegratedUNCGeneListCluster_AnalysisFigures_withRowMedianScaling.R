## 4/16/18
## Amy Olex
## Script that median centers the integrated PDX and TCGA data 
## before and after Limma normalization, then plots heatmaps and outputs TreeView data files. 
## UPDATE 12.15.21: To include new PDX samples.
## UPDATE 04.13.22: To include new PDX samples.

## Files I need: 
## Annotations and PAM50 classifications for the 151 PDX samples
## TCGA expression data: /Volumes/GoogleDrive/My Drive/000 OFFLINE LAPTOP FOLDERS/CCTR/Data/ChuckHarrell_BrainMetastasis_12-2016/2018.3.19_Integrating_PDX_with_TCGA_originFiles/TCGA_TPMvalues_fromTCGABioLinksQuery.RDS
## PAM50 gene list

source("http://www.bioconductor.org/biocLite.R")
library("NMF")
library("limma")
library("ggplot2")
library("ggdendro")
library("dendextend")
library("dplyr")

#setwd("~/Desktop/CCTR/Data/ChuckHarrell_BrainMetastasis_12-2016/2018.3.19_Integrating_PDX_with_TCGA/UNCExpandedSet/")
setwd("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq/2018-2022_merged_data")

#read in unnormalized
tcga <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2021_BulkRNASeq/2018-2021_merged_data/TCGA_Merge/Unnormalized_Log2TPM_UNCListGenes_TCGA817_and_Harrell_NOTHumanPercentCorrected_10MillionHumanReadFiltered.txt", row.names = 1, header=TRUE)
tcga <- tcga[,102:918]

#tcga_annot <- read.delim("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2021_BulkRNASeq/2018-2021_merged_data/TCGA_Merge/TCGA_subtypes.txt")
#row.names(tcga_annot) <- make.names(row.names(tcga_annot))

#tcga_filt <- tcga[,row.names(tcga_annot)]

## Row median scale data
unorm_rowMedians <- rowMedians(as.matrix(unnorm))
unnorm_scaled <- unnorm - unorm_rowMedians

limmanorm_rowMedians <- rowMedians(as.matrix(limmanorm))
limmanorm_scaled <- limmanorm - limmanorm_rowMedians


## Plot Heatmaps
metadata <- read.delim("../MetaData_TCGA-HarrellPDX-noMouse.txt")
row.names(metadata) <- metadata$Sample.ID

unnorm_scaled <- unnorm_scaled[,make.names(row.names(metadata)), drop=FALSE]
limmanorm_scaled <- limmanorm_scaled[,make.names(row.names(metadata)), drop=FALSE]

all(make.names(row.names(metadata))==names(limmanorm_scaled))
all(names(unnorm_scaled)==make.names(row.names(metadata)))
all(names(unnorm_scaled)==names(limmanorm_scaled))


### Get clusters and proportion counts
annot_data <- data.frame(PAM50=as.character(metadata$PAM50), Batch=as.character(metadata$Batch))
#my_colors <- lapply(lapply(annot_data, function(x){levels(as.factor(x))}), function(x){if(length(x)<=6){seq(from=1, to=length(x))}else{rainbow(length(x))}     })

#> levels(metadata$PAM50)
#[1] "Basal"  "Her2"   "LumA"   "LumB"   "Normal"

my_colors = list(PAM50=c("red", "magenta", "blue","cyan", "green3"), Batch=c("black", "white"))


###############
#### Get Unnormalized Clustering
###############

###### GENE EXPRESSION Unnormalized
jpeg(filename=paste("GeneExpHeatmap_UnNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_041818.jpeg", sep=""), width=20000, height=40000, res=2500, pointsize=5)
tree <- aheatmap(unnorm_scaled, distfun="correlation", hclustfun="ward", main=paste("UNCList Gene Expression\nUnNormalized Log2(TPM)\nRow Median Scaled BEFORE Clustering", sep=""), scale="none", annCol=annot_data, annColors=my_colors, 
                 fontsize=10, cexRow=4, treeheight=150, color=colorRampPalette(c('blue', 'white', 'red'))(20))
dev.off()

###### DENDROGRAM Unnormalized
clusters <- cutree(as.hclust(tree$Colv), k=4)
metadata$UnNormClusters <- clusters
#my_colors2 <- data.frame(PAM50=as.numeric(factor(metadata$PAM50)), Batch=as.numeric(factor(metadata$Batch)), Cluster=as.numeric(factor(metadata$UnNormClusters)))

my_colors2 <- recode(metadata$PAM50, "Basal"=2,"Her2"=6,"LumA"=4,"LumB"=5,"Normal"=3)


jpeg(filename=paste("ColumnDendrogram_UnNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_041818.jpeg", sep=""), width=10000, height=15000, res=1000, pointsize=1)

dend2 <- as.dendrogram(tree$Colv)
par(mai=c(0,0,1,4))
d1 <- color_branches(dend2, k=4, col = c("magenta","yellow","orange", "cyan"))
d1 %>% set("leaves_pch",15) %>%  # node point type
  set("leaves_cex", 1.5) %>%  # node point size
  set("leaves_col", my_colors2, order_value=TRUE) %>% #node point color
  plot(horiz = TRUE)
title(main="UNC List Gene Expr Sample Dendrogram\nNOT Normalized\nLog2(TPM) Row Median Centered", cex.main=20)
#legend(x=1500,y=900, cex = 10, bty = "n", title=expression(bold("Branch Colors")), legend=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4"), fill = c("magenta","cyan","yellow","orange") )
legend(x=2500,y=900, cex = 10, bty = "n", title=expression(bold("PAM50 Subtype (node colors)")), legend=c("Basal", "Her2","LumA","LumB","Normal"), fill = c("red", "pink", "darkblue","lightblue", "green") )

dev.off()


##### CORRELATION Matrix Unnormalized
cor_unnorm_scaled <- cor(unnorm_scaled)

## Have to do the below twice and overwrite the image to get the clusters correct
jpeg(filename=paste("CorHeatmap_UnNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_041818.jpeg", sep=""), width=20000, height=20000, res=2000, pointsize=5)
tree3 <- aheatmap(cor_unnorm_scaled, distfun="euclidean", hclustfun="ward", main=paste("UNC List Gene Expr Correlation\nUnNormalized Log2(TPM)\nRow Median Scaled Before Clustering", sep=""), scale="none", annCol=annot_data, annRow=annot_data, annColors=my_colors, 
                 fontsize=10, cexRow=4, treeheight=150, color=colorRampPalette(c('blue', 'white', 'red'))(20))
dev.off()


clusters <- cutree(as.hclust(tree3$Colv), k=4)
annot_data <- data.frame(PAM50=as.character(metadata$PAM50), Batch=as.character(metadata$Batch), Cluster=as.character(clusters))
my_colors$Cluster = c("magenta","cyan","yellow","orange")

jpeg(filename=paste("CorHeatmap_UnNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_withClusters_041818.jpeg", sep=""), width=10000, height=10000, res=1000, pointsize=1)
aheatmap(cor_unnorm_scaled, distfun="euclidean", hclustfun="ward", main=paste("UNC List Gene Expr Correlation\nUnNormalized Log2(TPM)\nRow Median Scaled", sep=""), scale="none", annCol=annot_data, annRow=annot_data, annColors=my_colors, 
         fontsize=10, cexRow=4, treeheight=150, color=colorRampPalette(c('blue', 'white', 'red'))(20))
dev.off()


### Convert to TreeView format
#hc <- as.hclust(tree$Colv)
#hc$dist.method <- "pearson"
#hc$method <- "wardD2"
#hr <- as.hclust(tree$Rowv)
#hr$dist.method <- "pearson"
#hr$method <- "wardD2"
#r2atr(hc, file="Unnormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianCentered.atr", distance=hc$dist.method, dec='.', digits=5)
#r2gtr(hr, file="Unnormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianCentered.gtr", distance=hr$dist.method, dec='.', digits=5)
#r2cdt(hr,hc,unnorm_scaled,labels=FALSE, description=FALSE, file="Unnormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianCentered.cdt", dec='.')




###########
#### Get Limma Normalized Clustering
###########

##### GENE EXPRESSION Limma Normalized

jpeg(filename=paste("GeneExpHeatmap_LimmaNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_041818.jpeg", sep=""), width=20000, height=40000, res=2500, pointsize=5)
tree4 <- aheatmap(limmanorm_scaled, distfun="correlation", hclustfun="ward", main=paste("UNCList Gene Expression\nLimmaNormalized Log2(TPM)\nRow Median Scaled BEFORE Clustering", sep=""), scale="none", annCol=annot_data, annColors=my_colors, 
                 fontsize=10, cexRow=4, treeheight=150, color=colorRampPalette(c('blue', 'white', 'red'))(20))
dev.off()

##### DENDROGRAM Limma Normalized

clusters <- cutree(as.hclust(tree4$Colv), k=4)
metadata$LimmaNormClusters <- clusters

jpeg(filename=paste("ColumnDendrogram_LimmaNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_041818.jpeg", sep=""), width=10000, height=15000, res=1000, pointsize=1)

dend2 <- as.dendrogram(tree4$Colv)
par(mai=c(0,0,1,4))
d1 <- color_branches(dend2, k=4, col = c("magenta","yellow","orange", "cyan"))
d1 %>% set("leaves_pch",15) %>%  # node point type
  set("leaves_cex", 1.5) %>%  # node point size
  set("leaves_col", my_colors2, order_value=TRUE) %>% #node point color
  plot(horiz = TRUE)
title(main="UNC List Gene Expr Sample Dendrogram\nLimma Normalized\nLog2(TPM) Row Median Centered", cex.main=20)
#legend(x=1500,y=900, cex = 10, bty = "n", title=expression(bold("Branch Colors")), legend=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4"), fill = c("magenta","cyan","yellow","orange") )
legend(x=2500,y=900, cex = 10, bty = "n", title=expression(bold("PAM50 Subtype (node colors)")), legend=c("Basal", "Her2","LumA","LumB","Normal"), fill = c("red", "pink", "darkblue","lightblue", "green") )

dev.off()


##### CORRELATION Matrix Limma Normalized
cor_limmanorm_scaled <- cor(limmanorm_scaled)
annot_data2 <- data.frame(PAM50=as.character(metadata$PAM50), Batch=as.character(metadata$Batch))

## Have to do the below twice and overwrite the image to get the clusters correct
jpeg(filename=paste("CorHeatmap_LimmaNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_041818.jpeg", sep=""), width=20000, height=20000, res=2000, pointsize=5)
tree5 <- aheatmap(cor_limmanorm_scaled, distfun="euclidean", hclustfun="ward", main=paste("UNC List Gene Expr Correlation\nLimmaNormalized Log2(TPM)\nRow Median Scaled Before Clustering", sep=""), scale="none", annCol=annot_data2, annRow=annot_data2, annColors=my_colors, 
                  fontsize=10, cexRow=4, treeheight=150, color=colorRampPalette(c('blue', 'white', 'red'))(20))
dev.off()


clusters2 <- cutree(as.hclust(tree5$Colv), k=4)
annot_data2 <- data.frame(PAM50=as.character(metadata$PAM50), Batch=as.character(metadata$Batch), Cluster=as.character(clusters2))
my_colors$Cluster = c("magenta","cyan","yellow","orange")

jpeg(filename=paste("CorHeatmap_LimmaNormalized_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianScaledBefore_withClusters_041818.jpeg", sep=""), width=10000, height=10000, res=1000, pointsize=1)
aheatmap(cor_limmanorm_scaled, distfun="euclidean", hclustfun="ward", main=paste("UNC List Gene Expr Correlation\nLimmaNormalized Log2(TPM)\nRow Median Scaled", sep=""), scale="none", annCol=annot_data2, annRow=annot_data2, annColors=my_colors, 
         fontsize=10, cexRow=4, treeheight=150, color=colorRampPalette(c('blue', 'white', 'red'))(20))
dev.off()



### Convert to TreeView format
#hc2 <- as.hclust(tree2$Colv)
#hc2$dist.method <- "pearson"
#hc2$method <- "wardD2"
#hr2 <- as.hclust(tree2$Rowv)
#hr2$dist.method <- "pearson"
#hr2$method <- "wardD2"
#r2atr(hc2, file="LimmaNorm_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianCentered.atr", distance=hc$dist.method, dec='.', digits=5)
#r2gtr(hr2, file="LimmaNorm_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianCentered.gtr", distance=hr$dist.method, dec='.', digits=5)
#r2cdt(hr2,hc2,limmanorm_scaled,labels=FALSE, description=FALSE, file="LimmaNorm_HarrellPDX-vs-TCGA817_UNCListGenesOnly_RowMedianCentered.cdt", dec='.')


#write.table(limmanorm_scaled, file="LimmaNorm_Log2TPM_UNCListGenes_TCGA817_and_Harrell_RowMedianScaled.txt",sep="\t", quote=FALSE)
#write.table(unnorm_scaled, file="Unnormalized_Log2TPM_UNCListGenes_TCGA817_and_Harrell_RowMedianScaled.txt",sep="\t", quote=FALSE)






