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

if(!dir.exists(here("data","plots"))) {
  dir.create(here("data","plots"), recursive = TRUE)
}

if(!dir.exists(here("data","SCEs","merge"))) {
  dir.create(here("data","SCEs","merge"), recursive = TRUE)
}

if(!dir.exists(here("data","SCEs","integration"))) {
  dir.create(here("data","SCEs","integration"), recursive = TRUE)
}

celltype <- "RecM"

# SCTransform, merge, and PCA --------------------------------------------------
out <- here("data","SCEs","merge",paste0("bal_",celltype,".SCT.obj.list.RData"))
s <- here("data","SCEs","merge",paste0("bal_",celltype,".merged.SEU.qs"))

if(!file.exists(out)) {
  seu <- qread(here("data","SCEs","manual.annot",paste0("bal_",celltype,".SEU.qs")),nthreads=16)
  DefaultAssay(seu) <- "RNA"
  seu <- DietSeurat(seu, assays = "RNA", dimreducs = NULL, graphs = NULL)

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

} else {if (!file.exists(s)) {
  # read object list
  obj.list <- qread(obj, nthreads=16)

  # find most variable features across samples to integrate
  var.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures=3000)

  # merge normalized samples
  seu <- merge(obj.list[[1]], y=obj.list[2:length(obj.list)], merge.data=TRUE)
  DefaultAssay(seu) <- "SCT"

  # manually set variable features of merged Seurat object
  VariableFeatures(seu) <- var.features

  # PCA
  seu <- RunPCA(seu, verbose=TRUE)

  qsave(seu, file=s, nthreads=16)

} else {
  seu <- qread(s, nthreads=16)
}}

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
out <- here("data","SCEs","integration",paste0("bal_",celltype,".integrated.SEU.qs"))
if(!file.exists(out)) {
  qsave(seu, file = out, nthreads=16)
}

# Leiden clustering ------------------------------------------------------------
if(!dir.exists(here("data","SCEs","clustering"))) {
  dir.create(here("data","SCEs","clustering"), recursive = TRUE)
}
out <- here("data","SCEs","clustering",paste0("bal_",celltype,".integrated.clustered.SEU.qs"))

if(!file.exists(out)) {

  seu <- FindNeighbors(seu, reduction="harmony",dims=1:30)

  plan("multicore", workers=20)
  print("Running Leiden clustering...")
  seu <- FindClusters(seu, algorithm = 4,
                      method = "igraph",
                      resolution = seq(0.2,3,0.2))

  # run UMAP
  seu <- RunUMAP(seu, assay="SCT", reduction = "harmony", dims=1:30)

  # save clustered object
  qsave(seu, file = out, nthreads=16)

} else {
  print("Using integrated object to perform re-clustering...")
  seu <- qread(here("data","SCEs","integration",paste0("bal_",celltype,".integrated.SEU.qs")), nthreads=16)
  seu <- FindNeighbors(seu, reduction = "harmony", dims=1:50)

  plan("multicore", workers=24)
  print("Running Leiden clustering...")
  seu <- FindClusters(seu, algorithm = 4,
                      method = "igraph",
                      resolution = seq(0.2,3,0.2))

  # run UMAP
  seu <- RunUMAP(seu, assay="SCT", reduction = "harmony", dims=1:50)

  # PrepSCTFindMarker
  plan("multicore", workers=4)
  seu <- PrepSCTFindMarkers(seu)

  # save clustered object
  qsave(seu,
        file = here("data","SCEs","clustering",paste0("bal_",celltype,".integrated.clustered.SEU.qs")),
        nthreads=16)

}

