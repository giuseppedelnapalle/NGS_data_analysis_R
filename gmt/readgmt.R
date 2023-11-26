#!/usr/bin/env Rscript
# source https://github.com/jhrcook/readgmt/blob/master/R/readgmt.R

#' Read a GMT file from MSigDB
#'
#' @description Reads a GMT file containing a gene set from
#'     \href{http://software.broadinstitute.org/gsea/index.jsp}{MSigDB} and
#'     returns a list of named gene sets
#'
#' @param file path to the downloaded GMT file
#' @param tidy turn into a tidy tibble with \code{tidy_gmt}
#'
#' @return returns a named list of the gene set
#'
#' @examples
#' kegg_path <- system.file("extdata", "c2.cp.kegg.v6.2.symbols.gmt.txt",
#'                          package = "readgmt")
#' kegg <- read_gmt(kegg_path)
#' head(names(kegg))
#' head(kegg$KEGG_COLORECTAL_CANCER)
#'
#' @export read_gmt
read_gmt <- function(file, tidy = FALSE) {
  con <- file(file, "r")
  gmt_lines <- readLines(file, warn = FALSE)
  close(con)
  rlist <- purrr::map(gmt_lines, parse_gmt_lines)
  rlist_names <- purrr::map_chr(gmt_lines, get_gmt_names)
  names(rlist) <- rlist_names
  if (tidy) rlist <- tidy_gmt(rlist)
  return(rlist)
}


#' Get genes of a gene set
#'
#' @description Returns a list of the genes from a single line in the GMT file
#'
#' @param gl a line of the GMT file
#'
#' @return a list of the genes
#'
#' @export parse_gmt_lines
parse_gmt_lines <- function(gl) {
  gl <- unlist(stringr::str_split(gl, "\\\t"), recursive = FALSE)
  gl <- gl[3:length(gl)]
  return(gl)
}


#' Get gene set name
#'
#' @description Returns the gene set name from a line of the GMT file
#'
#' @param gl a line of the GMT file
#'
#' @return the first value in the list (i.e. the gene set name)
#'
#' @export get_gmt_names
get_gmt_names <- function(gl) {
  gl <- unlist(stringr::str_split(gl, "\\\t"))
  gl_name <- gl[[1]]
  return(gl_name)
}
