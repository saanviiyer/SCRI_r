# load libraries
library(Seurat)
library(dplyr)

# read RDS
data <- readRDS("/Users/saanviiyer/Downloads/seurat_object_E9.5.rds")

# percent mitochondria count
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# violin plot to see quality control
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filtering bad data
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalization
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# top 10 genes that account for the highest variability
top10 <- head(VariableFeatures(data), 10)

# generating plots
plot1 <- VariableFeaturePlot(data)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,raster=TRUE)

all.genes <- rownames(data) # rows=genes, cols=cells
data <- SCTransform(data, verbose = FALSE)