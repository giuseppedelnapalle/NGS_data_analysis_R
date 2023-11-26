#!/usr/bin/env Rscript
# $title identify & remove outlying samples based on clustering
# $input gene expression data
# $author giuseppe
# $date Apr 2020

library(WGCNA)
library(flashClust)

# remove outliers with the Euclidean distance based sample network
remove_outliers <- function(mat = mat,  # genes as rows, samples as columns
                           phenoDT = phenoDT,
                           trait = trait,  # a vector
                           thresh_z = -2.5,
                           height = 10, width = 16, point_size = 14,
                           directory = directory,
                           file_suffix = ""){
  
  traitDt <- binarizeCategoricalColumns.forPlots(phenoDT[,trait], 
                                                 convertColumns = trait)
  
  A=adjacency(mat,type="distance")
  # this calculates the whole network connectivity
  k=as.numeric(apply(A,2,sum))-1
  # standardized connectivity
  Z.k=scale(k)
  
  # Designate samples as outlying
  # if their Z.k value is below the threshold
  # the color vector indicates outlyingness (red)
  outlierColor=ifelse(Z.k<thresh_z,"red","black")
  
  # calculate the cluster tree using flahsClust or hclust
  sampleTree = flashClust(as.dist(1-A), method = "average")
  
  # Convert traits to a color representation: 
  # where red indicates high values 
  traitColors=data.frame(numbers2colors(traitDt,signed=FALSE)) 
  dimnames(traitColors)[[2]]=paste(colnames(traitDt),"C",sep="") 
  datColors=data.frame(outlierC=outlierColor,traitColors) 
  
  # Plot the sample dendrogram and the colors underneath.
  pdf(file = paste0(directory,"Sample dendrogram and trait heatmap_",file_suffix,".pdf"),
      width = width, height = height, pointsize = point_size)
  plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                      colors=datColors,main="Sample dendrogram and trait heatmap")
  dev.off()
  
  # Remove outlying samples from expression and trait data 
  remove.samples= Z.k<thresh_z | is.na(Z.k) 
  print(paste("# removed samples",sum(remove.samples)))
  
  mat2 <- mat[,!remove.samples]
  return(mat2)
}
