#!/usr/bin/env Rscript

# This script performs HLCA v2 annotation in R v.4.2.2

# load packages ----------------------------------------------------------------
suppressPackageStartupMessages({
library(Seurat)
library(Azimuth)
})
set.seed(1990)

out <- "/group/canc2/anson/working/cf-eti-bal/data/SCEs/annotation/bal.annotated.SEU.rds"

if(!file.exists(out)) {
  dir.create("/group/canc2/anson/working/cf-eti-bal/data/SCEs/annotation", recursive=TRUE)
  # read object ------------------------------------------------------------------
  seu <- readRDS("bal_combined.preprocessed.SEU.rds")
  
  # run Azimuth ------------------------------------------------------------------
  seu <- RunAzimuth(seu,
                    "/group/canc2/anson/reference/annotation/HLCA_v2",
                    umap.name = 'azimuth.umap',
                    verbose = TRUE,
                    assay = "RNA",
                    k.weight = 50,
                    n.trees = 20,
                    mapping.score.k = 100)
  
  # write output - ---------------------------------------------------------------
  # save object
  saveRDS(seu, out)
}

print("HLCA v2 annotation Done!")
