#!/usr/bin/env Rscript
# functions for data frame operation

#' inner join multiple data frames (observations as columns and features as rows)
#' first column is for feature names
inner_join_multi <- function(df_lst, by_column, all_x_y=FALSE) {
  if (length(df_lst) < 2) {
    stop("At least two data frames are required for an inner join.")
  }
  
  df_nm <- names(df_lst)
  is_null <- is.null(df_nm)
  
  res <- df_lst[[1]]
  if (!is_null) {
    colnames(res)[2:ncol(res)] <- paste(df_nm[1], colnames(res)[2:ncol(res)], sep = "-")
    }
  
  for (i in 2:length(df_lst)) {
    new_df <- df_lst[[i]]
    if (!is_null) {
      colnames(new_df)[2:ncol(new_df)] <- paste(df_nm[i], colnames(new_df)[2:ncol(new_df)], sep = "-")
    }
    res <- merge(res, new_df, by = by_column, all = all_x_y)
  }
  
  return(res)
}
