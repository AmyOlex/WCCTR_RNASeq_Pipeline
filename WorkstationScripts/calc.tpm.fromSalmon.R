#' Calculates the TPM values for RNA-seq data from the summarized data output by Salmon and tximport.
#'
#' This function uses the method described at () to calculate the TPM (Transcripts Per Million) values from RNA-seq read counts and effective transcript lengths output by Salmon.
#'
#' @param data The list of matricies produced by tximport() of Salmon data.  Should contain a column named "length" and a column named "counts"
#' @param sample_names The list of names for each column.  Can pass in the list of files that was input to tximport().
#' @keywords tpm, gene, expression
#' @return A data frame of TPM expression values.
#' @examples
#' tpm <- calc.tpm.fromSalmon(data, colnames)
#' @export calc.tpm.fromSalmon
#' @author Amy L. Olex \email{alolex@@vcu.edu}
#'
calc.tpm.fromSalmon <- function(data, sample_names){

    RPK <- matrix(0, nrow=dim(data$counts)[1], ncol=dim(data$counts)[2])
    
    for(row in 1:dim(data$counts)[1]){
      for(col in 1:dim(data$counts)[2]){
        RPK[row,col] <- data$counts[row,col]/data$length[row,col]
      }
    }
    
    ##Calculate the sums of each column and divide by 1000000
    scale_factor <- colSums(RPK)/1000000
    
    ##Now divide all values in each column by the scaling factor
    TPM <- t(t(RPK)/scale_factor)
    colnames(TPM) <- sample_names
    row.names(TPM) <- row.names(data$counts)
    return(as.data.frame(TPM))
}
