#install.packages("Seurat")
library(Seurat)

#install.packages("dplyr")
library(dplyr)

#install.packages("presto")
#install.packages("devtools")
library(devtools)

basepath <- "/Users/saanviiyer/Downloads/filtered_gene_bc_matrices 2/hg19/"

pbmc.data <- Read10X(data.dir = basepath)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# to index: pbmc@meta.data[1:5, 1:2]

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)


plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,raster=TRUE)

plot1 + plot2

# scaling
all.genes <- rownames(pbmc) # rows=genes, cols=cells
pbmc <- ScaleData(pbmc, features = all.genes)

# runPCA option1
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimHeatmap(pbmc, dims = 1:8, cells = 500, balanced = TRUE)

# elbow
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:8)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.2)
pbmc <- FindClusters(pbmc, resolution = 1.2)

pbmc <- RunUMAP(pbmc, dims = 1:8)
# DimPlot(pbmc, reduction = "umap")
# DimPlot(pbmc, group.by = "RNA_snn_res.0.5")

# pbmc <- RunPCA(pbmc, dims = 1:8)
DimPlot(pbmc, reduction = "pca")
DimPlot(pbmc, reduction = "pca", group.by = "RNA_snn_res.1.2")

saveRDS(pbmc, paste0(basepath, "pbmc_seuratobj.RDS"))

# differential gene expression
# - explaining the difference between two gene populations through gene expression (turned on or off)

Idents(pbmc) <- "RNA_snn_res.0.2"

pbmc.markers <- FindAllMarkers(pbmc)

pbmc.markers %>%
  group_by(cluster) %>% # grouping dataframe by cluster column
  dplyr::filter(avg_log2FC > 1) %>% # filtering all groups by log2FC > 1
  slice_head(n = 10) %>% # selecting the top 10 genes for each cluster
  ungroup() -> top10 # restacking each group (AKA cluster) into one variable dataframe (i.e. top10)

# heatmap of the top10 upregulated genes per cluster res 0.2
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# gene ontology analysis
# ontology = way of creating categories of a specific area

# feature plot
naive_cd4_tcell_markers <- c("IL7R", "CCR7")
FeaturePlot(pbmc, features = naive_cd4_tcell_markers, reduction = "umap")
DimPlot(pbmc, group.by = "RNA_snn_res.0.5")

naive_cd4_markers_list <- list("naive_cd4_tcell_markers" = naive_cd4_tcell_markers)

# adding a module score based on the naive cd4 gene markers
pbmc <- AddModuleScore(pbmc, features = list("naive_cd4_tcell_markers" = naive_cd4_tcell_markers), name = "naive_cd4_tcell")
colnames(pbmc@meta.data)[9] <- "naive_cd4_tcell"
FeaturePlot(pbmc, features = "naive_cd4_tcell")
