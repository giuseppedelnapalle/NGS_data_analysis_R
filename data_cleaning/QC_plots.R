#!/usr/bin/env Rscript
# $title plots for sample-level QC
# $description
#  functions to create plots to assess sample-level data quality
# $author giuseppe
# $created Apr 2020

# load libraries
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(reshape2)

# heatmap of sample correlations
cor_heatmap <- function(mat = mat, # genes as rows, samples as columns
                        height = 8, width = 8, point_size = 14,
                        directory = ".", # no "/" trailing the path
                        name_cor = "", # a string
                        file_suffix = ""){
  cor_mat <- cor(mat)
  col <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
  
  # file name
  fn <- paste0(paste("correlation_heatmap", name_cor, file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  
  pdf(file = fn, width = width, height = height, pointsize = point_size)
  image(cor_mat, col = col, zlim = c(-1,1))
  dev.off()
}

# heatmap of the sample-to-sample distances
heatmap_sample_hclust <- function(mat = mat, # genes as rows, samples as columns
                                  annotation_row = NA, annotation_col = NA, # a data frame
                                  sample_col = NA, # a string
                                  show_rownames = T, show_colnames = T,
                                  labels_row = NULL, labels_col = NULL,
                                  height = 8.27, width = 11.69, # A4 landscape in inches
                                  directory = ".", file_suffix = ""){
  # calculate sample-to-sample distances
  dst <- dist(t(mat), method = "euclidean")
  dst_mat <- as.matrix(dst)
  color <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  # if a data frame is passed to annotation_row or annotation_col, a string specifying sample col should be given
  annot_called <- is.data.frame(annotation_row) | is.data.frame(annotation_col)
  stop_or_not <- ifelse(annot_called, !is.character(sample_col), F)
  if (stop_or_not) stop("sample_col not provided")
  
  fn <- paste0(paste("heatmap_hierarchical_clustering", file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  # fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  
  if (is.data.frame(annotation_row)) {
    r_nm <- annotation_row[, sample_col]
    # rename rows of annotation_row as sample_col if not
    if (sum(rownames(annotation_row) != r_nm) == nrow(annotation_row)) {
      rownames(annotation_row) <- r_nm
    }
    c_nm <- colnames(annotation_row)
    # suppress annotation with sample names
    annotation_row <- annotation_row[, !colnames(annotation_row) == sample_col]
    # convert annotation_row back to a data frame if it is coerced to a character vector
    if (!is.data.frame(annotation_row)) {
      annotation_row <- as.data.frame(annotation_row)
      rownames(annotation_row) <- r_nm
      colnames(annotation_row) <- c_nm[!c_nm == sample_col]
    }
  }
  
  if (is.data.frame(annotation_col)) {
    r_nm <- annotation_col[, sample_col]
    # rename rows of annotation_col as sample_col if not
    if (sum(rownames(annotation_col) != r_nm) == nrow(annotation_col)) {
      rownames(annotation_col) <- r_nm
    }
    c_nm <- colnames(annotation_col)
    # suppress annotation with sample names
    annotation_col <- annotation_col[, !colnames(annotation_col) == sample_col]
    if (!is.data.frame(annotation_col)) {
      annotation_col <- as.data.frame(annotation_col)
      rownames(annotation_col) <- r_nm
      colnames(annotation_col) <- c_nm[!c_nm == sample_col]
    }
  }
  
  pheatmap(mat = dst_mat,
           clustering_distance_rows = dst, clustering_distance_cols = dst,
           color = color,
           annotation_row = annotation_row, annotation_col = annotation_col,
           show_rownames = show_rownames, show_colnames = show_colnames,
           filename = fn, width = width, height = height,
           labels_row = labels_row, labels_col = labels_col)
  # annotation_row
  # The rows in the data and in the annotation are matched using corresponding row names
}

# density plot of gene expression values
density_plot <- function(mat = mat, # genes as rows, samples as columns. # samples <= 20
                         x_axis = "gene expression", legend_label = "sample",
                         alpha = .1,
                         height = 8, width = 18,
                         directory = ".", file_suffix = ""){
  # check # samples
  message(paste0("number of samples is ", ncol(mat)))
  if (ncol(mat) > 20) warning("for clarity, number of samples should be less than or equal to 20")
  
  # reshape data for ggplot 
  # append gene_id col
  mat$gene_id <- rownames(mat)
  # mat <- mat[,c(ncol(mat), 1:(ncol(mat)-1))]
  # melting
  md <- melt(mat, id="gene_id") # require reshape2 package
  colnames(md)[c(2, 3)] <- c(legend_label, x_axis)
  
  fn <- paste0(paste("density_plot", file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  # fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  
  plot <- ggplot(data = md, aes_(x = as.name(x_axis), fill = as.name(legend_label), color = as.name(legend_label))) +
    geom_density(alpha = alpha)
  ggsave(filename = fn, plot = plot, height = height, width = width, units = "cm")
}

# box plot of gene expression values (or other profiles) by sample
boxplot_by_sample <- function(mat = mat, # genes as rows, samples as columns
                              meta = meta, # a data frame
                              x = x, # a string specifying which col of meta to be used as x-axis
                              y = "gene expression",
                              fill = "cornflowerblue", # a string, either a colour name or a column name of meta
                              alpha = .5,
                              mask_x = TRUE, # logical, if TRUE then x-axis is labeled as "sample"
                              height = 21.01, width = 29.69, # A4 landscape
                              directory = ".", file_suffix = ""){
  # check # samples
  message(paste0("number of samples is ", ncol(mat)))
  if (ncol(mat) > 100) warning("for clarity, number of samples should be less than or equal to 100")
  
  # check if sample names in mat and meta match
  colnames(mat) %in% meta[, x]
  if (sum(colnames(mat) %in% meta[, x]) != ncol(mat)) stop("at least one sample name in mat not found in meta")
  
  fn <- paste0(paste("box_plot", file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  # fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  
  # set colour inside the boxes in two ways, depending on the fill argument
  if (fill %in% colnames(meta)) {
    # reshape data for ggplot
    meta <- meta[, c(x, fill)] # remove unnecessary columns
    rownames(meta) <- meta[, x]
    dt <- merge(meta, t(mat), by = "row.names") # inner join
    dt <- dt[, -1]
    md <- melt(dt, id = c(x, fill)) # require reshape2 package
    
    colnames(md)[ncol(md)] <- y
    if (mask_x) {
      x <- "sample"
      colnames(md)[1] <- x
    }
    
    plot <- ggplot(md, aes_(x = as.name(x), y = as.name(y), fill = as.name(fill))) +
      geom_boxplot(position = position_dodge(width = .1, preserve = "single"), color = "black", alpha = alpha)
    # alternative aes_string(x = x, y = y, fill = fill)
    ggsave(filename = fn, plot = plot, height = height, width = width, units = "cm")
  } else {
    meta <- as.data.frame(meta[, x])
    colnames(meta) <- x
    rownames(meta) <- meta[, x]
    dt <- merge(meta, t(mat), by = "row.names")
    dt <- dt[, -1]
    md <- melt(dt, id = x)
    
    colnames(md)[ncol(md)] <- y
    if (mask_x) {
      x <- "sample"
      colnames(md)[1] <- x
    }
    
    plot <- ggplot(md, aes_(x = as.name(x), y = as.name(y))) +
      geom_boxplot(position = position_dodge(width = .1, preserve = "single"), 
                   fill = fill, color = "black", alpha = alpha)
    ggsave(filename = fn, plot = plot, height = height, width = width, units = "cm")
  } # end of else
}

# high-density scatter plot
density_scatterplot <- function(mat = mat, # genes as rows, samples as columns
                                x = x,
                                y = y,
                                nrpoints = 100,
                                height = 8, width = 8, point_size = 14,
                                directory = ".", file_suffix = file_suffix){
  data_plt <- mat[,c(x, y)]
  xvar <- data_plt[,x]
  yvar <- data_plt[,y]
  
  # correlation coefficient
  r <- cor(xvar, yvar)
  
  # correlation test
  test <- cor.test(xvar, yvar)
  p <- test$p.value
  
  fn <- paste0(paste("scatterplot_smoothed_density", y, x, file_suffix, sep = "_"), ".pdf")
  # delete redundant "_" in file name
  # fn <- gsub("__*", "_", fn, fixed = F)
  fn <- gsub("_.pdf$", ".pdf", fn, fixed = F)
  fn <- paste(directory, fn, sep = "/")
  
  pdf(file = fn, width = width, height = height, pointsize = point_size)
  smoothScatter(data_plt, xlab = x, ylab = y, nrpoints = nrpoints,
                main = paste0(file_suffix, ", ", "r=", round(r, 3), ", ", "p=", p))
  dev.off()
}
