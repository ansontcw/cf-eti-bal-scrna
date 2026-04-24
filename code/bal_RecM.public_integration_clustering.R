#!/usr/bin/env Rscript

# This script performs integration and clustering analysis

# Load packages
suppressPackageStartupMessages({
library(harmony)
library(Seurat)
library(here)
library(qs)
library(ggplot2)
library(clustree)
library(cowplot)
library(forcats)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pals)
})
set.seed(1990)
options(future.globals.maxSize = 500000000000, future.seed=TRUE)
date <- Sys.Date()

celltype <- "public"

if(!dir.exists(here("data","public_datasets","plots"))) {
  dir.create(here("data","public_datasets","plots"), recursive = TRUE)
}

# SCTransform, merge, and PCA --------------------------------------------------
out <- here("data","public_datasets","SCEs",paste0(celltype,".SCT.obj.list.RData"))
s <- here("data","public_datasets","SCEs",paste0(celltype,".merged.SEU.qs"))

if(!file.exists(out)) {
  seu <- qread(here("data","public_datasets","SCEs",paste0(celltype,".processed.SEU.qs")),nthreads=16)
  DefaultAssay(seu) <- "RNA"

  # split object by experiment
  obj.list <- SplitObject(seu, split.by="batchID")

  # normalization
  obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- SCTransform(x, vst.flavor = "v2", verbose=TRUE)
  })

  # save object list
  qsave(obj.list, file = out, nthreads=16)

  # find most variable features across samples to integrate
  var.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures=3000)

  # merge normalized samples
  seu <- merge(obj.list[[1]], y=obj.list[2:length(obj.list)], merge.data=TRUE)
  DefaultAssay(seu) <- "SCT"

  # manually set variable features of merged Seurat object
  VariableFeatures(seu) <- var.features

  # PCA
  seu <- RunPCA(seu, verbose=TRUE)

  # save merged object
  qsave(seu, file=s, nthreads=16)

  obj.list <- NULL

} else {
  seu <- qread(s, nthreads=16)
}

# Harmony integration ----------------------------------------------------------

# run Harmony
print("Running Harmony integration...")
DefaultAssay(seu) <- "SCT"
seu <- RunHarmony(seu,
                  group.by.vars=c("batchID"),
                  reduction="pca",
                  assay.use="SCT",
                  reduction.save="harmony",
                  plot_convergence=TRUE,
                  early_stop=FALSE,
                  max_iter=500)

# save integrated object
out <- here("data","public_datasets","SCEs",paste0(celltype,".integrated.SEU.qs"))
if(!file.exists(out)) {
  qsave(seu, file = out, nthreads=16)
}

# Leiden clustering ------------------------------------------------------------
if(!dir.exists(here("data","SCEs","clustering"))) {
  dir.create(here("data","SCEs","clustering"), recursive = TRUE)
}
out <- here("data","public_datasets","SCEs",paste0(celltype,".integrated.clustered.SEU.qs"))

if(!file.exists(out)) {
  seu <- FindNeighbors(seu, reduction="harmony",dims=1:50)

  plan("multicore", workers=20)
  print("Running Leiden clustering...")
  seu <- FindClusters(seu, algorithm = 4,
                      method = "igraph",
                      resolution = seq(0.2,3,0.2))

  # run UMAP
  seu <- RunUMAP(seu, assay="SCT", reduction = "harmony", dims=1:50)
  
  # PrepSCTFindmarkers
  plan("multicore", workers=6)
  seu <- PrepSCTFindMarkers(seu)
  
  # save clustered object
  qsave(seu, file = out, nthreads=16)
  
  # run clustree
  ggsave(plot=clustree(seu, prefix="SCT_snn_res."), device="pdf",
         width = 10,
         height = 30,
         file=here("data","public_datasets","plots",celltype,paste0(celltype,".clustree.",date,".pdf")))
}
