## Amy Olex
## 10.28.21
## Analysis to create initial heatmaps.  

library(ComplexHeatmap)
library(RNASeqBits)
library(ggplot2)
library(ggrepel)

setwd("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2022_Bulk_RNASeq/")

## Import Quantile Normalized Log2 TPM Data with zero exp genes removed
human_TPM <- read.delim("processed_data/log2TPM_human_ZeroExpRemoved_QuantileNorm_2022_BulkRNASeq_04.13.22.tsv", header=TRUE, row.names=1, stringsAsFactor=FALSE)
mouse_TPM <- read.delim("processed_data/log2TPM_mouse_ZeroExpRemoved_QuantileNorm_2022_BulkRNASeq_04.13.22.tsv", header=TRUE, row.names=1, stringsAsFactor=FALSE)

## import sample metadata
samples <- read.delim("sample_metadata_master_04.13.22.csv", header=TRUE, stringsAsFactors = FALSE, sep=",")
row.names(samples) <- make.names(samples$Sample.ID)


human_TPM_genes <- human_TPM[,"symbol",drop=FALSE]
mouse_TPM_genes <- mouse_TPM[,"symbol",drop=FALSE]

human_TPM.norm <- human_TPM[,row.names(samples)]
mouse_TPM.norm <- mouse_TPM[,row.names(samples)]

#######
## normalize
### I already imported the quantile normalized values.
########
## Upper Quantile Normalize the imported Log2 TPM values
## this uses all data in the quantile function, including those that are not expressed
#data.quantileAll <- apply(human_TPM, 2, function(x){quantile(x, 0.75)})
#human_TPM.norm <- as.data.frame(t(t(human_TPM) / data.quantileAll))

#data.quantileAll_mouse <- apply(mouse_TPM, 2, function(x){quantile(x, 0.75)})
#mouse_TPM.norm <- as.data.frame(t(t(mouse_TPM) / data.quantileAll_mouse))


##### PCA using all samples ####
## PCA Analysis of selected PDX samples for PAM50 genes
human_TPM.norm.scaled <- as.data.frame(t(scale(t(as.matrix(human_TPM.norm)))))

pca <-  prcomp(t(human_TPM.norm.scaled))

## The percent variance Explained:
data.frame(summary(pca)$importance)[, 1:min(5, ncol(summary(pca)$importance))]


  
png(filename = "figs/PCAplot_log2TPM_human_ZeroExpRemoved_QuantileNorm_2022_BulkRNASeq_04.13.22.png", res = 150, width=1200, height=1000)  
  
  colorby <- "PDX.line"
  ggplot(data = data.frame(pca$x, samples, samples = samples$PDX.line, stringsAsFactors = F), 
             aes(x = as.numeric(PC1), y = as.numeric(PC2), label = Sample.ID)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca)$importance[2,2] * 100, digits = 2), "% variability" ))

dev.off()


####### ExpressionHeatmap ######

centered_colors <- center.palette(human_TPM.norm.scaled, palette_length = 100)

gene_vars = apply(human_TPM.norm.scaled, 1, var)

top2000 = human_TPM.norm.scaled[names(sort(gene_vars, decreasing = TRUE))[1:2000],,drop=FALSE]

gene_sym_top2000 <- human_TPM_genes[row.names(top2000),,drop=FALSE]
all(row.names(gene_sym_top2000)==row.names(top2000))
row.names(top2000) <- gene_sym_top2000$symbol

centered_colors <- center.palette(top2000, palette_length = 100)
#my_annot <- HeatmapAnnotation(df = annot_data_50m, col = list(Tissue=c("Brain"="red", "Liver" = "blue", "Lung" = "black") ) )

png(filename=paste("figs/Heatmap_top2000variable_log2TPM_human_ZeroExpRemoved_QuantileNorm_2022_BulkRNASeq_04.13.22.png", sep=""), width=2000, height=3000, res=150)

Heatmap(as.matrix(top2000), col = centered_colors$colors,  
        column_dend_height = unit(2, "cm"), show_row_dend = TRUE, row_dend_width = unit(2, "cm"),
        clustering_method_columns="ward.D2", clustering_method_rows="ward.D2",
        show_row_names=FALSE, show_column_names=TRUE) 

dev.off()


########################################### Not using below code yet

## Row Median Center Genes before plotting
## Row median scale data
#rowMedians <- rowMedians(as.matrix(human_TPM_pam50_selected))
#human_TPM_pam50_selected_scaled <- human_TPM_pam50_selected - rowMedians

#centered_colors <- center.palette(human_TPM_pam50_selected_scaled, palette_length = 100)

#jpeg(filename=paste("2018.07.20_Heatmap_Human_SelectedPDX_PAM50.jpeg", sep=""), width=1500, height=1300, res=150)
#  aheatmap(t(human_TPM_pam50_selected_scaled), distfun="correlation", hclustfun="ward", main=paste("Selected PDX MGT Samples\nPAM50 Gene Set Human Expression\nUpper-Quantile Norm, Row-Median Centered", sep=""), scale="none", treeheight=100, color = centered_colors)
#dev.off()


## PCA Analysis of selected PDX samples for PAM50 genes
pca <-  prcomp(t(scale(as.matrix(human_TPM_pam50_selected))))

## The percent variance Explained:
data.frame(summary(pca)$importance)[, 1:min(5, ncol(summary(pca)$importance))]

colorby <- "PAM50.Subtype" # covariates[2]
pt <- ggplot(data = data.frame(pca$x, samples_selected, samples = samples_selected$PDX.line, stringsAsFactors = F), 
             aes(x = as.numeric(PC1), y = as.numeric(PC2), label = PDX.line)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA with batch, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 4.5) +
  scale_color_manual(values=c("red", "magenta", "cyan")) +
  geom_text_repel(colour = "black", size = 3, point.padding = unit(0.25, "lines")) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  labs(color = colorby) +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca)$importance[2,2] * 100, digits = 2), "% variability" )) +
  theme_bw()

png(filename=paste("Figure_6B_2018.07.25_PCA_Plot_Human_SelectedPDX_PAM50.png", sep=""), width=2000, height=1200, res=250)
plot(pt)
dev.off()


########
### Figure 6C
##################### >50% species Correlation Heatmaps ######################
#### Only using samples with >50% mouse

mouse_greater50_samples <- row.names(samples)[samples$X.Mouse > .50]

## Note the above are not exactly mutually exclusive because there are 4 samples with less than 50% in both human and mouse due to unmapped read percentage being high.

mouse_TPM_greater50 <- mouse_TPM.norm[,mouse_greater50_samples]

samples_mouse_greater_50 <- samples[mouse_greater50_samples,]

#### Plot Mouse
#annot_data_50m <- data.frame(Tissue=samples_mouse_greater_50$tissue, Batch=as.character(samples_mouse_greater_50$batch), PDX.line=samples_mouse_greater_50$PDX.line)
annot_data_50m <- data.frame(Tissue=factor(samples_mouse_greater_50$Tissue))
annot_data_50m[annot_data_50m$Tissue == "NA Liver",] <-  "Liver"
annot_data_50m$Tissue <- factor(annot_data_50m$Tissue)
my_colors_50m <- lapply(lapply(annot_data_50m, function(x){levels(as.factor(x))}), function(x){if(length(x)<=6){seq(from=1, to=length(x))}else{rainbow(length(x))}     })

#annot_data_50m$human.percent <- samples_mouse_greater_50$human.percent
#annot_data_50m$mouse.percent <- samples_mouse_greater_50$mouse.percent


cor_mouse_greater_50 <- cor(mouse_TPM_greater50)

#jpeg(filename=paste("2018.07.20_Heatmap_mouse_GreaterThan50PercentMouse_allGenes.jpeg", sep=""), width=2500, height=2300, res=150, pointsize=5)
#aheatmap(cor_mouse_greater_50, distfun="euclidean", hclustfun="average", scale="none", annCol=annot_data_50m,annRow=annot_data_50m, annColors=my_colors_50m, 
#         fontsize=10, treeheight=200, color = c("white","red"))
#dev.off()




#png(filename=paste("Figure_6C_2018.07.20_GeneExpressionHeatmap_mouse_GreaterThan50PercentMouse_Top2000VariableGenes_withWHIM30.png", sep=""), width=2000, height=1200, res=200, pointsize=5)

#aheatmap(top2000, distfun="correlation", hclustfun="ward.D2", scale="none", annCol=annot_data_50m,annColors=my_colors_50m, 
#         fontsize=10, cexRow=0, cexCol=0, treeheight=50, color = centered_colors)
#dev.off()






## SKIPPING!
####################
###### Figure 6B: PAM50 selected MGT PDX samples Gene Expression Heatmap
####################

## Import Gene Subsets
pam50 <- read.delim("~/Desktop/CCTR/Data/ChuckHarrell_BrainMetastasis_12-2016/8-18-17_DataFreeze_InitialHeatmaps/PAM50_geneSym_ensembl.txt", header=FALSE, stringsAsFactors = FALSE)
names(pam50) <- "symbol"

human_TPM_pam50 <- human_TPM.norm[which(human_TPM_genes$symbol %in% pam50$symbol),]

PAM50_genes <- human_TPM_genes[which(human_TPM_genes$symbol %in% pam50$symbol),,drop=FALSE]

## Put PAM50 genes in order and replace the row names with gene symbols
human_TPM_pam50 <- human_TPM_pam50[row.names(PAM50_genes),,drop = FALSE]

## sanity check
all(row.names(human_TPM_pam50) == row.names(PAM50_genes))

## replace gene IDs with gene symbols
row.names(human_TPM_pam50) <- PAM50_genes$symbol



