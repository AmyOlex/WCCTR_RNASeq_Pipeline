#' Calculates the TPM values for RNA-seq data from the summarized data output by FeatureCounts
#'
#' This function uses the method described at () to calculate the TPM (Transcripts Per Million) values from RNA-seq read counts and effective transcript lengths output by FeatureCounts.
#'
#' @param metadata The dataframe with the meta data for each gene, must include a gene length column.
#' @param data The dataframe with ONLY counts in it and the row.names are the same as the metadata table.
#' @keywords tpm, gene, expression
#' @return A data frame of TPM expression values.
#' @examples
#' tpm <- calc.tpm.fromFeatureCounts(metadata, data)
#' @export calc.tpm.fromFeatureCounts
#' @author Amy L. Olex \email{alolex@@vcu.edu}
#'
calc.tpm.fromFeatureCounts <- function(metadata, data){
  
  RPK <- matrix(0, nrow=dim(data)[1], ncol=dim(data)[2])
  
  for(row in 1:dim(data)[1]){
    for(col in 1:dim(data)[2]){
      RPK[row,col] <- data[row,col]/metadata$Length[row]
    }
  }
  
  ##Calculate the sums of each column and divide by 1000000
  scale_factor <- colSums(RPK)/1000000
  
  ##Now divide all values in each column by the scaling factor
  TPM <- t(t(RPK)/scale_factor)
  colnames(TPM) <- colnames(data)
  row.names(TPM) <- rownames(data)
  return(as.data.frame(TPM))
}
