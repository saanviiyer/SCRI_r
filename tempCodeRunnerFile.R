# ORGANOGENESIS DATASET---------------------------------------------------------

# load libraries----------------------------------------------------------------
library(Seurat)
library(dplyr)
library(devtools)
library(ggplot2)

basepath <- "/Users/saanviiyer/Downloads/org_data/"

#create counts matrix-----------------------------------------------------------
cts <- ReadMtx(mtx = "/Users/saanviiyer/Downloads/org_data/matrix.mtx.gz",
               features = "/Users/saanviiyer/Downloads/org_data/features.tsv.gz",
               cells = "/Users/saanviiyer/Downloads/org_data/barcodes.tsv.gz")
org_data <- CreateSeuratObject(counts = cts, project = "pmc3k", min.cells = 3, min.features = 200)

# percent mitochondria count----------------------------------------------------
org_data[["percent.mt"]] <- PercentageFeatureSet(org_data, pattern = "^MT-")

# visualize QC metrics as a violin plot-----------------------------------------
VlnPlot(org_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
org_data <- subset(org_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize---------------------------------------------------------------------
org_data <- NormalizeData(org_data, normalization.method = "LogNormalize", scale.factor = 10000)

org_data <- FindVariableFeatures(org_data, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(org_data), 10)
plot1 <- VariableFeaturePlot(org_data)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# scaling-----------------------------------------------------------------------
all.genes <- rownames(org_data) # rows=genes, cols=cells
org_data <- ScaleData(org_data, features = all.genes)

# run PCA-----------------------------------------------------------------------
org_data <- RunPCA(org_data, features = VariableFeatures(object = org_data))

DimHeatmap(org_data, dims = 1:11, cells = 500, balanced = TRUE)

ElbowPlot(org_data)

# produce separate resolutions--------------------------------------------------
org_data <- FindNeighbors(org_data, dims = 1:11)
org_data <- FindClusters(org_data, resolution = 0.5)
org_data <- FindClusters(org_data, resolution = 0.2)
org_data <- FindClusters(org_data, resolution = 1.2)

org_data <- RunUMAP(org_data, dims = 1:11)

# plot dimensionality-----------------------------------------------------------
DimPlot(org_data, reduction = "umap")
DimPlot(org_data, group.by = "RNA_snn_res.1.2")

#saveRDS(org_data, paste0(basepath, "organogenesis_seuratobj.RDS"))

# differential gene expression--------------------------------------------------
Idents(org_data) <- "RNA_snn_res.1.2"

org_data.markers <- FindAllMarkers(org_data)

org_data.markers %>%
  group_by(cluster) %>% # grouping dataframe by cluster column
  dplyr::filter(avg_log2FC > 1) %>% # filtering all groups by log2FC > 1
  slice_head(n = 10) %>% # selecting the top 10 genes for each cluster
  ungroup() -> top10 # restacking each group (AKA cluster) into one variable dataframe (i.e. top10)

# heatmap of the top10 upregulated genes per cluster res 0.2
DoHeatmap(org_data, features = top10$gene) + NoLegend() 


# gene ontology analysis--------------------------------------------------------

#  feature plot
markers <- c("MEOX1", "FABP5")
FeaturePlot(org_data, features = markers, reduction = "umap")
# -> this produces an error: Error in UseMethod(generic = "FetchData", object = object) : no applicable method for 'FetchData' applied to an object of class "data.frame"
# debugger says that error is about whether org_data is a Seurat Object
DimPlot(org_data, group.by = "RNA_snn_res.0.5")

markers_list <- list("markers" = markers)

# adding a module score based on the gene markers
org_data <- AddModuleScore(org_data, features = list("markers" = markers), name = "cell")
colnames(org_data@meta.data)[9] <- "cell"
FeaturePlot(org_data, features = "cell", min.cutoff = "q1") + ggtitle(paste(markers[1],"and",markers[2]))

FeaturePlot(org_data, features = top10$gene[1:6])


