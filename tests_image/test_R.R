#!/opt/conda/bin/Rscript

# test loading all packages
library(MAST)
library(monocle)
library(scater)
library(scran)
library(Seurat)
library(SingleCellExperiment)

# print versions
sessionInfo()