## Amy Olex
## 01.30.25
## Script to check if the files and directories in a samplesheet exists or if there are errors.  
##
## This script checks to ensure all files pointed to in a sample sheet actualy exist or if there are errors.
##
## TODO: Need to expand this to be able to input any sample sheet with any column names and using tab or csv format so it is more versatile.

library(readr)


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Required. The path and name of the sample sheet file to validate.", metavar="character"),
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("A sample sheet must be provided.", call.=FALSE)
}

samplesheet <- opt$file

#samplesheet <- "/lustre/home/harrell_lab/bulkRNASeq/jaxpdxnf/sampleSheets/Harrell_Bulk_RNASeq_Master_nfJAX_SampleSheet_121924.csv"

data <- read.delim(samplesheet, header=TRUE, sep=",")

fq1 <- data[which(!file.exists(data$fastq_1)),]

fq2 <- data[which(!file.exists(data$fastq_2)),]

if(dim(fq1)[1] > 0 | dim(fq2)[1] > 0){
  message("Error: Some files do not exists.")
  print("Fastq1 Files:")
  print(fq1)
  print("Fastq2 Files:")
  print(fq2)
} else {
  print("All files validated!")
}

sessionInfo()