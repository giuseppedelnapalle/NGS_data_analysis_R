#!/usr/bin/env Rscript
# preprocessing functions for scRNA-seq data from H5 files

library(foreach)
library(doParallel)

# select cells based on cell identities
# data an m-by-n sparse Matrix of class "dgCMatrix", where m is # genes, n is # cells
# meta a factor
# select a character vector
select_cells <- function(data, meta, select){
  mat <- data[,colnames(data) %in% names(meta)]
  meta <- as.character(meta)
  mat_f <- mat[,meta %in% select]
}

# find cell identities for two objects
# mat an m-by-n sparse Matrix of class "dgCMatrix", where m is # genes, n is # cells
find_cell_ident <- function(mat, meta, select,
                            mat_2, meta_2, select_2){
  c_n <- names(meta)[as.character(meta) %in% select]
  c_i <- as.character(meta[c_n])
  c_n_2 <- names(meta_2)[as.character(meta_2) %in% select_2]
  c_i_2 <- as.character(meta_2[c_n_2])
  c_i_m <- as.numeric(as.factor(c(c_i, c_i_2)))
}

# calculate number of features
calc_n_features <- function(mat){
  n_f <- apply(mat, 2, function(x) sum(x>0))
}

# calculate number of counts
calc_n_counts <- function(mat){
  n_c <- apply(mat, 2, sum)
}

# filter cells by number of features
filter_cells_n_f <- function(mat, thresh){
  n_f <- calc_n_features(mat)
  mat_f <- mat[,n_f >= thresh]
}

# filter cells by library size
filter_cells_lib_s <- function(mat, thresh){
  n_c <- calc_n_counts(mat)
  mat_f <- mat[,n_c >= thresh]
}

# filter genes
filter_genes <- function(mat, min_cells=1){
  f <- apply(mat, 1, function(x) sum(x>0) >= min_cells)
  mat_f <- mat[f,]
}

# match genes of two gene expression matrices
match_genes <- function(mat, mat_2){
  c_genes <- intersect(rownames(mat), rownames(mat_2))
  mat_f <- mat[rownames(mat) %in% c_genes,]
  mat_2_f <- mat_2[rownames(mat_2) %in% c_genes,]
  res <- list(mat_f, mat_2_f)
}

# a wrapper function of filter_genes and match_genes
filter_match_genes <- function(mat, mat_2, min_cells=1){
  mat_f <- filter_genes(mat, min_cells=min_cells)
  mat_2_f <- filter_genes(mat_2, min_cells=min_cells)
  lst <- match_genes(mat_f, mat_2_f)
}

# calculate min library size
# mat (m-by-n) m genes, n cells
min_lib_size <- function(mat){
  min(apply(mat, 2, sum))
}

# down sample gene counts
# new counts = expected counts given the adjusted lib size
# mat (m-by-n) m genes, n cells
down_sample_exp <- function(mat, adj_lib_s){
  mat_d <- apply(mat, 2, function(x){
    floor((x / sum(x)) * adj_lib_s)
    # (x / sum(x)) * adj_lib_s
  })
}

# # down sample gene counts
# # mat (m-by-n) m genes, n cells
# down_sample <- function(mat, adj_lib_s, iter=100, seed=64548){
#   n_g <- nrow(mat)
#   n_c <- ncol(mat)
#   set.seed(seed)
#   # mat_d <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#   mat_d <- sapply(1:n_c, function(i){
#     print(paste0("i: ", i))
#     c_it <- sapply(1:iter, function(j){
#       print(paste0("j: ", j))
#       g_c <- rep(0, n_g)
#       tb_smp <- table(sample.int(n_g, adj_lib_s, replace = TRUE, prob = mat[,i]/sum(mat[,i])))
#       v_smp <- as.vector(tb_smp)
#       g_c[match(as.integer(names(tb_smp)), c(1:n_g))] <- v_smp
#       return(g_c)
#     })
#     c_rd <- apply(c_it, 1, median)
#     return(c_rd)
#   })
#   return(mat_d)
# }

# down sample gene counts
# mat (m-by-n) m genes, n cells
# down_sample <- function(mat, adj_lib_s, iter=100, seed=64548){
#   n_g <- nrow(mat)
#   n_c <- ncol(mat)
#   set.seed(seed)
# 
#   # sapply method
#   # mat_d <- sapply(1:n_c, function(i){
#   #   print(paste0("i: ", i))
#   #   c_it <- sapply(1:iter, function(j){
#   #     print(paste0("j: ", j))
#   #     g_c <- rep(0, n_g)
#   #     tb_smp <- table(sample.int(n_g, adj_lib_s, replace = TRUE, prob = mat[,i]/sum(mat[,i])))
#   #     v_smp <- as.vector(tb_smp)
#   #     g_c[match(as.integer(names(tb_smp)), c(1:n_g))] <- v_smp
#   #     return(g_c)
#   #   })
#   #   c_rd <- apply(c_it, 1, median)
#   #   return(c_rd)
#   # })
# 
#   # for loop method
#   mat_d <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#   for (i in 1:n_c) {
#     print(paste0("i: ", i))
#     c_it <- matrix(0, n_g, iter)
#     for (j in 1:iter) {
#       print(paste0("j: ", j))
#       g_c <- rep(0, n_g)
#       tb_smp <- table(sample.int(n_g, adj_lib_s, replace = TRUE, prob = mat[,i]/sum(mat[,i])))
#       # v_smp <- as.vector(tb_smp)
#       g_c[match(as.integer(names(tb_smp)), c(1:n_g))] <- as.vector(tb_smp)
#       c_it[,j] <- g_c
#     }
#     mat_d[,i] <- apply(c_it, 1, function(x) round(median(x)))
#   }
#   return(mat_d)
# 
#   # mat_lst <- lapply(1:iter, function(i){
#   #   mat_d <- apply(mat, 2, function(x){
#   #     g_c <- rep(0, n_g)
#   #     tb_smp <- table(sample.int(n_g, adj_lib_s, replace = TRUE, prob = x/sum(x)))
#   #     v_smp <- as.vector(tb_smp)
#   #     g_c[match(as.integer(names(tb_smp)), c(1:n_g))] <- v_smp
#   #     return(g_c)
#   #   })
#   # })
#   # m_agg <- do.call(cbind, mat_lst)
#   # mat_d <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#   # for (j in 1:n_c) {
#   #   print(paste0("j: ", j))
#   #   # c_it <- m_agg[,seq(1, ncol(m_agg), n_c)]
#   #   mat_d[,j] <- apply(m_agg[,seq(1, ncol(m_agg), n_c)], 1, function(x) round(median(x)))
#   # }
#   #
#   # print(dim(mat_d))
#   # return(mat_d)
# 
# }

# library(parallel)
# down_sample <- function(mat, adj_lib_s, iter=100, seed=64548){
#   n_g <- nrow(mat)
#   n_c <- ncol(mat)
#   set.seed(seed)
#   
#   mat_lst <- mclapply(1:iter, function(i){
#     apply(mat, 2, function(x){
#       g_c <- rep(0, n_g)
#       tb_smp <- table(sample.int(n_g, adj_lib_s, replace = TRUE, prob = x/sum(x)))
#       v_smp <- as.vector(tb_smp)
#       g_c[match(as.integer(names(tb_smp)), c(1:n_g))] <- v_smp
#       return(g_c)
#     })
#   }, mc.cores = detectCores())
#   
#   m_agg <- do.call(cbind, mat_lst)
#   mat_d <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#   
#   for (j in 1:n_c) {
#     c_it <- m_agg[,seq(j, ncol(m_agg), n_c)]
#     mat_d[,j] <- apply(c_it, 1, function(x) round(median(x)))
#   }
#   
#   return(mat_d)
# }

# library(foreach)
# library(doParallel)

# down sample gene counts
# mat (m-by-n) m genes, n cells
down_sample <- function(mat, adj_lib_s, iter = 100, seed = 64548, num_cores = detectCores() - 1) {
  n_g <- nrow(mat)
  n_c <- ncol(mat)
  set.seed(seed)
  
  cl <- makeCluster(num_cores) # initialize cluster
  registerDoParallel(cl) # register cluster for parallel processing
  
  # foreach method with rbind
  mat_lst <- foreach(j = 1:n_c) %dopar% {
    # print(paste0("j: ", j))
    c_it <- matrix(0, n_g, iter)
    for (k in 1:iter) {
      # print(paste0("k: ", k))
      g_c <- rep(0, n_g)
      tb_smp <- table(sample.int(n_g, adj_lib_s, replace = TRUE, prob = mat[, j] / sum(mat[, j])))
      g_c[match(as.integer(names(tb_smp)), c(1:n_g))] <- as.vector(tb_smp)
      c_it[, k] <- g_c
    }
    apply(c_it, 1, function(x) round(median(x)))
  }
  
  stopCluster(cl) # stop cluster
  
  mat_d <- do.call(cbind, mat_lst)
  rownames(mat_d) <- rownames(mat)
  colnames(mat_d) <- colnames(mat)
  return(mat_d)
}
