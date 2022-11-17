# Amy Olex
# 11/17/22
# Collecting metrics of analyzed single cell data.

library("optparse")

debug_file <- "~/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/SingleCell/debug_files/metrics_summary.csv"

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis. Be sure to include whether this is before or after filtering out mouse or dead cells.", metavar="character"),
  make_option(c("-c", "--configfile"), type="character", 
              help="Required. Input config file with the sample name in the first column and the path to the metric summary file in the second column.", metavar="character"),
  make_option(c("-x", "--cellranger_report"), type="logical", 
              help="If the report files are the metric_summary.csv file from 10X, then we don't need to enter in the file name (default = FALSE)", 
              default = FALSE, action = "store_true", metavar="logical"),
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

if (is.null(opt$configfile)){
  print_help(opt_parser)
  stop("A configuration file with sample information must be provided.", call.=FALSE)
}

runID <- opt$runid
inFile <- opt$configfile
reportDir <- opt$outdir
reportName <- paste0(reportDir, runID, "_MetricSummaryReport.txt")
cellranger <- opt$cellranger_report
  
  
#runID <- "DEBUG"
#inFile <- "./debug_files/debug_config.txt"
#reportDir <- "./debug_files/"
#reportName <- paste0(reportDir, runID, "_MetricSummaryReport.txt")

#### Ok, print out all the current running options as a summary.

print("Summary of input options:\n")
print(paste("Run Name:", runID))
print(paste("ConfigFile:", inFile))
print(paste("Report Output:", reportName))
print(paste("Processing 10X report:", cellranger))

# Read in the provided config file and loop for each row.
toProcess = read.table(inFile, header=FALSE, sep=",")
print(paste(dim(toProcess)[1], " rows were found."))

# Initialize summary report data frame
print(paste("Processing row 1 from sample", toProcess[1,1]))

if(cellranger){
  reportData <- read.csv(paste0(as.character(toProcess[1,2]), "/metrics_summary.csv"))
  report_sampleIDs <- as.character(toProcess[1,1])
}else{
  reportData <- read.csv(as.character(toProcess[1,2]))
  report_sampleIDs <- as.character(toProcess[1,1])
}



for(i in 2:dim(toProcess)[1]){
  print(paste("Processing row", i, "from sample", toProcess[i,1]))
  
  report_sampleIDs <- append(report_sampleIDs, as.character(toProcess[i,1]))
  
  if(cellranger){
    reportData <- rbind(reportData, read.csv(paste0(as.character(toProcess[i,2]), "/metrics_summary.csv")))
  }else{
    reportData = rbind(reportData, read.csv(as.character(toProcess[i,2])))
  }
  
  
}

row.names(reportData) <- report_sampleIDs

write.table(reportData, file = reportName, quote = FALSE, row.names=TRUE)


print("Completed!")





