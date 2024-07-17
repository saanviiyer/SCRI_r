# load libraries
library(Seurat)
library(dplyr)

data <- readRDS("/Users/saanviiyer/Downloads/cerebellar_embryogenesis_seuratobject.RDS")

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
plot1 <- VariableFeaturePlot(data)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,raster=TRUE)

all.genes <- rownames(data) #rows = genes, cols=cells

data <- SCTransform(data, verbose = FALSE)

data <- ScaleData(data, features = all.genes)

data <- RunPCA(data, features = VariableFeatures(object = data))

DimHeatmap(data, dims = 1:8, cells = 500, balanced = TRUE)

ElbowPlot(data)