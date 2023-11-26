#!/usr/bin/env Rscript
# $title convert ensembl_ID to gene_symbol
# $author giuseppe
# $date Dec 2022


# 0 setting up R session --------------------------------------------------

# turn off stringsAsFactors
options(stringsAsFactors = F)

# directory of gene information
dir_gene <- "~/Project_Data/R_data/BioMart"

# 1 load data -------------------------------------------------------------

load(file = paste(dir_gene, "gene_info.RData", sep = "/"))


# 2 convert ensembl_ID to gene_symbol -------------------------------------

ensembl_to_geneSym <- function(mat = mat,
                               gene_mapping = gene_info){
  g0 <- rownames(mat)
  
  # remove duplicated rows in gene_mapping
  gene_mapping <- gene_mapping[,c("Gene.stable.ID", "Gene.name")]
  gene_mapping <- unique(gene_mapping)
  
  # confirm all genes of mat in gene_mapping
  if (sum(g0 %in% gene_mapping$Gene.stable.ID) == nrow(mat)) {
    
    rn <- gene_mapping$Gene.name[match(g0, gene_mapping$Gene.stable.ID)]
    print(paste0("# unique gene_symbol ", length(unique(rn))))
    
    if (length(unique(rn)) == nrow(mat)) {
      rownames(mat) <- rn
      return(mat)
      
    } else {
      print("duplicated mapping found.")
      
      # select gene with max median across samples
      dup_gs <- unique(rn[duplicated(rn)])
      
      print("duplicated genes: ")
      print(dup_gs)
      print(paste0("# duplicated genes ", length(dup_gs)))
      
      keep_dup <- sapply(dup_gs, function(x){
        # index of duplicated genes
        dup_ens_x <- gene_mapping[gene_mapping$Gene.name == x,]$Gene.stable.ID
        dup_ens_x <- g0[g0 %in% dup_ens_x]
        dup_idx <- match(dup_ens_x, g0)
        medians <- apply(mat[dup_idx,], 1, median)
        # gene index with max median
        keep <- dup_idx[order(medians, decreasing = T)][1]
        return(keep)
      })

      # duplicated genes
      dup_ens <- gene_mapping[gene_mapping$Gene.name %in% dup_gs,]$Gene.stable.ID
      dup_ens <- g0[g0 %in% dup_ens]
  
      keep <- (! rownames(mat) %in% dup_ens) | c(1:nrow(mat)) %in% keep_dup
      mat2 <- mat[keep,]
      rownames(mat2) <- gene_mapping$Gene.name[match(rownames(mat2), gene_mapping$Gene.stable.ID)]
      return(mat2)
    }  # end of else (length(unique(rn)))
  } else {
    print("some genes of mat not in gene_mapping.")
  } # to be finished when necessary
}
