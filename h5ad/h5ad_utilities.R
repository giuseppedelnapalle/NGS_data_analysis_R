#!/usr/bin/env Rscript
# utilities to manage H5AD files

library(anndata)
library(Seurat)

#' obtain a dataset from an AnnData object
#' @param ann an AnnData object
#'        path_dataset path to the dataset in the object, e.g. "raw/X"
#' @return a sparse of class dgRMatrix
get_ds_anndata <- function(ann, path_dataset){
  p <- strsplit(path_dataset, "/")[[1]] # a character vector
  ds <- recursive_extract(ann, p)
}

#' define a function that recursively extracts nested elements of a list
recursive_extract <- function(x, s) {
  if (length(s) == 0) {
    return(x)
  } else {
    return(recursive_extract(x[[s[1]]], s[-1]))
  }
}

#' read cell annotation file
read_cell_annot <- function(filename){
  n <- max(count.fields(filename, sep = ','))
  tb <- read.table(filename, sep = ",", header = FALSE, 
                      col.names = paste0("V", seq_len(n)), fill = TRUE)
  cell_ann <- lapply(1:nrow(tb), function(i) {
    x <-unname(unlist(tb[i,]))
    x <- x[x != ""]
    x <- x[2:length(x)]
  })
  names(cell_ann) <- tb[,1]
  return(cell_ann) # a list of vectors
}

#' create filters to select cells of interest
#' @param various number of lists specifying how cells are filtered,
#' e.g. list("obs/cell_ontology_class", c("t cell", "b cell"))
#' @return a transformed list of vectors
create_filters <- function(...){
  args <- list(...)
  f <- lapply(1:length(args), function(i) unlist(args[[i]][2]))
  names(f) <- unlist(lapply(1:length(args), function(i) args[[i]][1]))
  return(f)
}

#' filter cells
filter_cells <- function(ann, filters){
  f <- lapply(1:length(filters), function(i){
    v <- as.vector(get_ds_ann(ann, names(filters)[i]))
    v %in% filters[[i]]
  })
  s <- bw_and(f)
  ann[s,]
}

bw_and <- function(lst){
  res <- lst[[1]]
  if (length(lst) > 1) {
    for (i in 2:length(lst)) {
      res <- res & lst[[i]]
    }
  }
  return(res)
}
