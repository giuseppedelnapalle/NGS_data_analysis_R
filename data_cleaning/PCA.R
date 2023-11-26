#!/usr/bin/env Rscript
# $title PCA
# $description
#  functions for PCA and visualisation
# $author giuseppe
# $date Apr 2020

# load libraries
library(ggplot2)
library(ggfortify)

# plot variance explained by PCs
plot_var_explained <- function(singular_value = singular_value, # a matrix, e.g. s$v (matrix whose columns contain the right singular vectors)
                               height = 8, width = 8, point_size = 14,
                               directory = ".",
                               file_suffix = ""){
  fn <- paste0(paste("variance_explained", file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  # fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  pdf(file = fn, width = width, height = height, pointsize = point_size)
  plot(singular_value^2/sum(singular_value^2),
       main = "Variance explained by PCs",
       xlab = "PC",
       ylab = "Variance_explained")
  dev.off()
}

# PCA scatter plot
# sample order of mat, meta, batch should match
# require ggfortify package
PCA_scatterplot <- function(mat = mat,  # samples as rows, genes as columns
                     meta = meta, # a data frame
                     batch = batch, # a data frame
                     color = color, # a string specifying col to be used for labeling samples w/ colors
                     scale = F,
                     loadings = F,
                     loadings_label = F,
                     loadings_label_size = 3,
                     height = 10, width = 16, point.size = 2, # point.size for size argument of geom_point
                     directory = ".",
                     file_suffix = ""){
  print("sample order of mat, meta, batch should match.")
  
  # compute principal components
  pca_res <- prcomp(mat, scale. = scale)
  
  df_col <- cbind(meta, batch)
  df_pca <- cbind(mat, as.data.frame(df_col[,color]))
  colnames(df_pca)[ncol(df_pca)] <- color
  
  plot = autoplot(pca_res, data = df_pca, colour=color, 
                  main = "PCA scatter plot",
                  loadings =  loadings, loadings.label = loadings_label, 
                  loadings.label.size = loadings_label_size, size = point.size)
  
  fn <- paste0(paste("PCA_scatter_plot", color, file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  # fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  ggsave(filename = fn, plot = plot, height = height, width = width, units = "cm")
}

# sample order of mat_singular_v, meta, batch should match
PCA_scatterplot2 <- function(mat_singular_v = mat_singular_v, # a matrix, e.g. s$v (matrix whose columns contain the right singular vectors)
                     components = components,    # a vector
                     meta = meta,
                     batch = batch,
                     color = color,
                     height = 10, width = 16, point.size = 2, # point.size for size argument of geom_point
                     directory = ".",
                     file_suffix = ""){
  print("sample order of mat_singular_v, meta, batch should match.")
  
  i <- components[1]
  j <- components[2]
  
  data_plt <- cbind(data.frame(PCx = mat_singular_v[,i], PCy = mat_singular_v[,j]), meta, batch)
  data_plt <- data_plt[,c("PCx", "PCy", color)]
  
  plot = ggplot(data_plt, aes_string(x = "PCx", y = "PCy", color = color)) +
    geom_point(size = point.size) + 
    labs(title="PCA scatter plot", x=paste0("PC", i), y=paste0("PC", j))
  
  fn <- paste0(paste("PCA_scatter_plot", color, file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  # fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  ggsave(filename = fn, plot = plot, height = height, width = width, units = "cm")
}
