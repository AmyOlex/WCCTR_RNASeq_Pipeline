# Amy Olex
# 11/18/22
# Collecting metrics of analyzed single cell data from the HTML reports.

library("optparse")
library("rvest")
library("xml2")

debug_file <- "~/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/SingleCell/debug_files/web_summary.html"

# Read and parse HTML file
doc.html = htmlTreeParse(debug_file, useInternal = TRUE)

doc.html = rvest::read_html(debug_file)

txt <- doc.html %>% html_nodes("body") %>% html_text()

grep("\"rows\":\\[\\[\"Number of Reads\",\"", txt)


