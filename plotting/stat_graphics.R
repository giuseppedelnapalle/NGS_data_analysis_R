#!/usr/bin/env Rscript
# statistical graphics

library(ggplot2)

#' create density plots for a list of vectors
density_plot <- function(vector_list, dir_fig, group_val,
                         x_lab="n_counts", group_lab="shape", fill=TRUE, 
                         alpha=.3, width=12, height=10, unit="cm", dpi=300){
  if (length(vector_list) != length(group_val)) {
    stop("length of vector_list is not equal to length of group_val")
  }
  
  grp_lst <- lapply(1:length(vector_list), function(i){
    rep(group_val[i], length(vector_list[[i]]))
  })
  df <- data.frame(unlist(vector_list), as.factor(unlist(grp_lst)))
  colnames(df) <- c(x_lab, group_lab)
  
  if (fill) {
    ggp <- ggplot(data = df, aes(x=.data[[x_lab]], fill=.data[[group_lab]]))
  } else {
    ggp <- ggplot(data = df, aes(x=.data[[x_lab]], color=.data[[group_lab]]))
  }
  p <- ggp + geom_density(alpha=alpha)
  
  fn <- paste0(paste("density_plot", x_lab, group_lab,  sep = "_"), ".pdf")
  ggsave(fn, plot = p, path = dir_fig, width = width, height = height, units = unit, dpi = dpi)
}
