#!/usr/bin/env Rscript

#' Transpose a dgRMatrix and simultaneously convert it to dgCMatrix
#' @param inmat input matrix in dgRMatrix format
#' @return A dgCMatrix that is the transposed dgRMatrix
#' @export transpose_dgRMatrix
transpose_dgRMatrix <- function(inmat) {
    if(class(inmat) != 'dgRMatrix')
        stop('inmat is not of class dgRMatrix')
    out <- new('dgCMatrix',
               i=inmat@j,
               p=inmat@p,
               x=inmat@x,
               Dim=rev(inmat@Dim),
               Dimnames=rev(inmat@Dimnames)
               )
    out
}
