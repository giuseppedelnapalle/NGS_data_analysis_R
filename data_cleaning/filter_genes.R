#!/usr/bin/env Rscript
# $title remove lowly expressed genes in gene expression data
# $input gene expression data
# $author giuseppe
# $date Apr 2020

# remove genes if TPM < thresh_tpm in over thresh_prop of samples
filter_genes <- function(mat = mat, # genes as rows, samples as columns
                         thresh_tpm = 1,
                         thresh_prop = .9){
  del <- rowSums(mat < thresh_tpm)/ncol(mat) > thresh_prop
  mat_flt <- mat[!del,]
  return(mat_flt) # return a matrix 
}
