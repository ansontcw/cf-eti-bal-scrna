#!/usr/bin/env Rscript

# This script performs GO BP enrichment analysis on selected cell types
# Run this script on directory: /group/canc2/anson/working/cf-et-bal

# Load packages
suppressPackageStartupMessages({
library(Seurat)
library(here)
library(qs)
library(forcats)
library(ggplot2)
library(rstatix)
library(dplyr)
library(DOSE)
library(tibble)
library(msigdbr)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(stringr)
library(readr)
})

set.seed(1990)
options(future.globals.maxSize = 500000000000, future.seed=TRUE)
date <- Sys.Date()

# CF vs HC ---------------------------------------------------------------------
# celltypes <- setNames(c("Epithelial", "RecM", "TRM-CCL"),c("Epithelial", "RecM", "TRM-CCL"))
# compare.list <- setNames(c("CFM_vs_HC","CFS_vs_HC"),c("CFM_vs_HC","CFS_vs_HC"))

# ETI
celltypes <- setNames(c("TRM-CCL","RecM"),
                      c("TRM-CCL","RecM"))

compare.list <- setNames(c("CF.ETI vs CF.ETI.BL",
                           "CF.UT vs CF.UT.BL",
                           "CF.ETI vs HC",
                           "CF.UT vs HC"),
                         c("CF.ETI_vs_CF.ETI.BL",
                           "CF.UT_vs_CF.UT.BL",
                           "CF.ETI_vs_HC",
                           "CF.UT_vs_HC"))
# Run GOBP enrichment ----------------------------------------------------------
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

bigdf <- data.frame()
for (celltype in names(celltypes)) {
  # loop through each comparison
  for (compare in names(compare.list)) {
    # if(dir.exists(here("data","GOBP",celltypes[celltype],compare, "GO_BP"))) {
    #   up.df <- read.csv(here("data","GOBP",celltypes[celltype],compare, "GO_BP",paste0(celltypes[celltype],".",compare,".GO_BP-up.csv")))
    #   down.df <- read.csv(here("data","GOBP",celltypes[celltype],compare, "GO_BP",paste0(celltypes[celltype],".",compare,".GO_BP-down.csv")))
    #   # append df to list
    #   bigdf <- rbind(bigdf,rbind(up.df,down.df))
    # }
    # # Create directory
    # else {
  if(!dir.exists(here("data", "GOBP", celltype, compare))) {
    dir.create(here("data", "GOBP", celltype, compare,"GO_BP"), recursive = TRUE)
  }
  
  # read DE files
  up <- read.csv(here("data","DE",celltypes[celltype],compare,"DEG",paste0(celltypes[celltype],".",compare,".DEG-up.csv")))
  down <- read.csv(here("data","DE",celltypes[celltype],compare,"DEG",paste0(celltypes[celltype],".",compare,".DEG-down.csv")))
  # get ensemblDB
  getCols <- setNames(c("SYMBOL","ENTREZID"),c("SYMBOL","ENTREZID"))
  
  genes <- data.frame(
    lapply(getCols, function(column) {
      mapIds(
        x = org.Hs.eg.db,
        keys = unique(c(up$gene,down$gene)),
        keytype = "SYMBOL",
        column = column)
    }),
    row.names = unique(c(up$gene,down$gene)))
  
  # GO pathway enrichment
  go.up.bp <- enrichGO(genes[up$gene,]$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                       pAdjustMethod = "BH", readable=TRUE, pvalueCutoff = 0.01)

  go.down.bp <- enrichGO(genes[down$gene,]$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                       pAdjustMethod = "BH", readable=TRUE, pvalueCutoff = 0.01)
  
  if (!is.null(go.up.bp)) {
    # Remove redundant GO terms
    go.up.bp <- clusterProfiler::simplify(go.up.bp)
    
    # write results
    go.up.bp@result %>%
      dplyr::rename(Pathway=Description) %>%
      dplyr::mutate(Celltype = names(celltypes[celltype])) %>%
      dplyr::mutate(Comparison=compare.list[compare]) %>%
      dplyr::mutate(Direction="Up") %>%
      dplyr::mutate(Genes=gsub("/","; ",go.up.bp@result$geneID)) %>%
      dplyr::select(Pathway, Celltype, Comparison, Direction, ID, Count, pvalue, p.adjust, Genes) %>%
      slice_head(n = 10) -> up.df
      write_csv(up.df, file = here("data","GOBP",celltype,compare, "GO_BP",
                            paste0(celltype,".",compare,".GO_BP-up.csv")))
    # append df to list
    bigdf <- rbind(bigdf,rbind(up.df))
    }
  
  if (!is.null(go.down.bp)) {
    # Remove redundant GO terms
    go.down.bp <- clusterProfiler::simplify(go.down.bp)
    
    # write results
    go.down.bp@result %>%
      dplyr::rename(Pathway=Description) %>%
      dplyr::mutate(Celltype = names(celltypes[celltype])) %>%
      dplyr::mutate(Comparison=compare.list[compare]) %>%
      dplyr::mutate(Direction="Down") %>%
      dplyr::mutate(Genes=gsub("/","; ",go.down.bp@result$geneID)) %>%
      dplyr::select(Pathway, Celltype, Comparison, Direction, ID, Count, pvalue, p.adjust, Genes) %>%
      slice_head(n = 10) -> down.df
      write_csv(down.df, file = here("data","GOBP",celltype,compare, "GO_BP",
                            paste0(celltype,".",compare,".GO_BP-down.csv")))
      # append df to list
      bigdf <- rbind(bigdf,rbind(down.df))
    }
    # }
  }
}

# write the big df
# write.csv(bigdf, file=here("data","GOBP","CF_vs_HC.bigdf.csv"), row.names = FALSE)
write.csv(bigdf, file=here("data","GOBP","ETI.bigdf.csv"), row.names = FALSE)
