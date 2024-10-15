
files <- read.delim("TCGA-BRCA_fileList.txt", header=FALSE)
f1 <- read.delim(files[1,1], skip = 1)
f1 <- f1[,c("gene_id", "gene_name", "tpm_unstranded"), drop=FALSE]
sname <- unlist(strsplit(basename(files[1,1]), split = '\\.'))[1]
names(f1) <- c("gene_id", "gene_name", sname)

for(fname in files[2:length(files$V1),1]){
	sname <- unlist(strsplit(basename(fname), split = '\\.'))[1]
	f <- read.delim(fname, skip = 1)
	f <- f[,c("tpm_unstranded"), drop=FALSE]
	names(f) <- sname
	f1 <- cbind(f1,f)
}

write.table(f1, file="TCGA-BRCA_STAR-Count_TPM.csv", quote=FALSE, row.names=FALSE, sep=",")
