#!/usr/bin/env Rscript
# $title functions for gene set functional annotation
# $desctiption 
#   GO or KEGG enrichment analysis over multiple gene sets
# $author giuseppe
# $date Apr 2020

library(clusterProfiler)

# enrich GO & generate dot plot or bar plot
enrichGO_multi <- function(gene_list = gene_list,
                           gene_label = gene_label,
                           drop_label = NULL,
                           GO_type = c("CC", "MF", "BP"),
                           OrgDb = "org.Hs.eg.db",
                           keyType = "SYMBOL", # "ENSEMBL", "ENTREZID"
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           background_gene = NULL,
                           qvalueCutoff  = 0.1,
                           plot_type = "dot.plot", # "dot.plot" or "bar.plot"
                           n_showCategory = 10,
                           width = 10,
                           height = 8,
                           point.size_plt = 14,
                           directory = directory) {
  clst.ls = levels(as.factor(gene_label)) %>% as.character()
  # drop label if specified
  if (!is.null(drop_label)) clst.ls = clst.ls[!clst.ls == drop_label]
  # enrich GO loop
  eGO.ls = lapply(clst.ls, function(clst) {
    gene_c = gene_list[gene_label == clst]
    res_GO = lapply(GO_type, function(GO) {
      # enrich GO
      eGO = enrichGO(gene          = gene_c,
                     OrgDb         = OrgDb,
                     keyType       = keyType,
                     ont           = GO,
                     pAdjustMethod = pAdjustMethod,
                     pvalueCutoff  = pvalueCutoff,
                     universe      = background_gene,
                     qvalueCutoff  = qvalueCutoff)
      # plot
      if (plot_type == "dot.plot") {
        filename = paste0(directory, "dot_plot_GO_enrichment_", GO, "_", clst,".pdf")
        pdf(file = filename, width = width, height = height, pointsize = point.size_plt)
        plot_tmp = dotplot(eGO, showCategory=n_showCategory)
        print(plot_tmp)
        dev.off()
      } else {
        filename = paste0(directory, "bar_plot_GO_enrichment_", GO, "_", clst,".pdf")
        pdf(file = filename, width = width, height = height, pointsize = point.size_plt)
        plot_tmp = barplot(eGO, showCategory=n_showCategory)
        print(plot_tmp)
        dev.off()
      } # end of else
      return(eGO)
    }) # end of lapply GO
    return(res_GO)
  }) # end of lapply clst
  return(eGO.ls)
} # end of function

# enrich KEGG pathways & generate dot plot or bar plot
# keytypes(org.Hs.eg.db)
enrichKEGG_multi <- function(gene_list = gene_list,
                           gene_label = gene_label,
                           drop_label = NULL,
                           organism = "hsa",
                           OrgDb="org.Hs.eg.db",
                           keyType = "kegg",
                           fromType="SYMBOL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           background_gene = NULL,
                           qvalueCutoff  = 0.1,
                           plot_type = "dot.plot", # "dot.plot" or "bar.plot"
                           n_showCategory = 10,
                           width = 10,
                           height = 8,
                           point.size_plt = 14,
                           directory = directory) {
  clst.ls = levels(as.factor(gene_label)) %>% as.character()
  # drop label if specified
  if (!is.null(drop_label)) clst.ls = clst.ls[!clst.ls == drop_label]
  # enrich KEGG loop
  eKEGG.ls = lapply(clst.ls, function(clst) {
    gene_c = gene_list[gene_label == clst]
    gene.entrez = bitr(gene_c, fromType=fromType, toType="ENTREZID", OrgDb=OrgDb)
    gene.entrez = unique(gene.entrez[,2])
    eKEGG = enrichKEGG(gene = gene.entrez,
                       organism = organism,
                       keyType = keyType,
                       pAdjustMethod = pAdjustMethod,
                       pvalueCutoff  = pvalueCutoff,
                       universe      = background_gene,
                       qvalueCutoff  = qvalueCutoff)
    # plot
    if (plot_type == "dot.plot") {
      filename = paste0(directory, "dot_plot_KEGG_enrichment_", clst,".pdf")
      pdf(file = filename, width = width, height = height, pointsize = point.size_plt)
      plot_tmp = dotplot(eKEGG, showCategory=n_showCategory)
      print(plot_tmp)
      dev.off()
    } else {
      filename = paste0(directory, "bar_plot_KEGG_enrichment_", clst,".pdf")
      pdf(file = filename, width = width, height = height, pointsize = point.size_plt)
      plot_tmp = barplot(eKEGG, showCategory=n_showCategory)
      print(plot_tmp)
      dev.off()
    } # end of else
    return(eKEGG)
  }) # end of lapply clst
  return(eKEGG.ls)
} # end of function


