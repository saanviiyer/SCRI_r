# ORGANOGENESIS DATASET---------------------------------------------------------

# load libraries----------------------------------------------------------------
library(Seurat)
library(dplyr)
library(devtools)
library(ggplot2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

# install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)