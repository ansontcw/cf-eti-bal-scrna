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
celltypes <- setNames(c("TRM-CCL","Epithelial","RecM"),
                      c("TRM-CCL","Epithelial","RecM"))
#compare.list <- setNames(c("CF.ETI.BL vs CF.UT.BL","CF.ETI.BL vs HC","CF.UT.BL vs HC","CF.ETI vs CF.ETI.BL",
#                           "CF.UT vs CF.UT.BL","CF.ETI vs HC","CF.UT vs HC","CF.ETI vs CF.UT"),
#                         c("CF04-06.BL_vs_CF01-03.BL","CF04-06.BL_vs_HC","CF01-03.BL_vs_HC","CF.ETI_vs_CF04-06.BL",
#                           "CF.UT_vs_CF01-03.BL","CF.ETI_vs_HC","CF.UT_vs_HC","CF.ETI_vs_CF.UT"))

compare1 <- setNames(c("CFN vs HC",
                       "CFB vs HC",
                       "CFB vs CFN"),
                     c("CFN_vs_HC",
                       "CFB_vs_HC",
                       "CFB_vs_CFN"))

compare2 <- setNames(c("CF.ETI vs CF.ETI.BL",
                       "CF.ETI vs HC",
                       "CF.UT vs CF.UT.BL",
                       "CF.UT vs HC"),
                     c("CF.ETI_vs_CF.ETI.BL",
                       "CF.ETI_vs_HC",
                       "CF.UT_vs_CF.UT.BL",
                       "CF.UT_vs_HC"))

# Run GOBP enrichment ----------------------------------------------------------
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

bigdf <- data.frame()
for (celltype in names(celltypes)) {
  if (celltype == "TRM-CCL") {
    compare.list <- c(compare1,compare2)
  } else if(celltype == "Epithelial") {
    compare.list <- compare1
  } else {
    compare.list <- compare2
  }
  # loop through each comparison
  for (compare in names(compare.list)) {
    up <- read.csv(here("data","DE",celltypes[celltype],compare,"DEG",paste0(celltypes[celltype],".",compare,".DEG-up.csv")))
    down <- read.csv(here("data","DE",celltypes[celltype],compare,"DEG",paste0(celltypes[celltype],".",compare,".DEG-down.csv")))
    
    up %>% 
      dplyr::mutate(Direction="Up") %>%
      dplyr::mutate(Celltype=celltype) %>%
      dplyr::mutate(Gene=up$gene) %>%
      dplyr::mutate(Comparison=compare.list[compare]) %>%
      dplyr::select(Gene, Celltype, Comparison, Direction, avg_log2FC, pct.1, pct.2, p_val_adj) -> up
    
    down %>% 
      dplyr::mutate(Direction="Down") %>%
      dplyr::mutate(Celltype=celltype) %>%
      dplyr::mutate(Gene=down$gene) %>%
      dplyr::mutate(Comparison=compare.list[compare]) %>%
      dplyr::select(Gene, Celltype, Comparison, Direction, avg_log2FC, pct.1, pct.2, p_val_adj) -> down
    
    # append df to list
    bigdf <- rbind(bigdf,rbind(up,down))
  }
}

# write the big df
# write.csv(bigdf, file=here("data","GOBP","CF_vs_HC.bigdf.csv"), row.names = FALSE)
write.csv(bigdf, file=here("data","DE","TableE9.csv"), row.names = FALSE)
