# load libraries-----------------
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)

# load data----------------------
#developmental days
seuobj89 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj89.RDS")
seuobj94 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj94.RDS")
seuobj110 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj110.RDS")
seuobj115 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj115.RDS")
seuobj125 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj125.RDS")

#tumors
seuobj_atrt <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_atrt.RDS")
seuobj_dipg <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_dipg.RDS")
seuobj_mb <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_mb.RDS")

# visualization for developmental days---------
DimPlot(seuobj89, reduction = "umap", group.by = "Main_cluster_name")
DimPlot(seuobj89, reduction = "umap", group.by = "sub_cluster_name")

DimPlot(seuobj94, reduction = "umap", group.by = "Main_cluster_name")
DimPlot(seuobj94, reduction = "umap", group.by = "sub_cluster_name")

DimPlot(seuobj110, reduction = "umap", group.by = "Main_cluster_name")
DimPlot(seuobj110, reduction = "umap", group.by = "sub_cluster_name")

DimPlot(seuobj115, reduction = "umap", group.by = "Main_cluster_name")
DimPlot(seuobj115, reduction = "umap", group.by = "sub_cluster_name")

DimPlot(seuobj125, reduction = "umap", group.by = "Main_cluster_name")
DimPlot(seuobj125, reduction = "umap", group.by = "sub_cluster_name")

#Pre Processing for Tumors=============

#atrt pre processing
seuobj_atrt[["percent.mt"]] <- PercentageFeatureSet(seuobj_atrt, pattern = "^MT-")
seuobj_atrt <- subset(seuobj_atrt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seuobj_atrt <- NormalizeData(seuobj_atrt, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_atrt <- FindVariableFeatures(seuobj_atrt, selection.method = "vst", nfeatures = 2000)
seuobj_atrt <- ScaleData(seuobj_atrt, features = rownames(seuobj_atrt))
seuobj_atrt <- RunPCA(seuobj_atrt, features = VariableFeatures(object = seuobj_atrt))
DimHeatmap(seuobj_atrt, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(seuobj_atrt)
seuobj_atrt <- FindNeighbors(seuobj_atrt, dims = 1:15)
seuobj_atrt <- FindClusters(seuobj_atrt, resolution = 0.8)
seuobj_atrt <- RunUMAP(seuobj_atrt, dims = 1:15)
DimPlot(seuobj_atrt, reduction = "umap", group.by = "seurat_clusters")


#dipg pre processing
seuobj_dipg[["percent.mt"]] <- PercentageFeatureSet(seuobj_dipg, pattern = "^MT-")
seuobj_dipg <- subset(seuobj_dipg, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seuobj_dipg <- NormalizeData(seuobj_dipg, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_dipg <- FindVariableFeatures(seuobj_dipg, selection.method = "vst", nfeatures = 2000)
seuobj_dipg <- ScaleData(seuobj_dipg, features = rownames(seuobj_dipg))
seuobj_dipg <- RunPCA(seuobj_dipg, features = VariableFeatures(object = seuobj_dipg))
DimHeatmap(seuobj_dipg, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(seuobj_dipg)
seuobj_dipg <- FindNeighbors(seuobj_dipg, dims = 1:15)
seuobj_dipg <- FindClusters(seuobj_dipg, resolution = 0.8)
seuobj_dipg <- RunUMAP(seuobj_dipg, dims = 1:15)
DimPlot(seuobj_dipg, reduction = "umap", group.by = "seurat_clusters")

#mb pre processing
DimHeatmap(seuobj_mb, dims = 1:18, cells = 500, balanced = TRUE)
DimPlot(seuobj_mb, reduction = "umap", group.by = "seurat_clusters")
ElbowPlot(seuobj_mb)


#Integration of tumors===========
#creating tumor seuobj list
tumors.seuobj.list <- list("atrt" = seuobj_atrt,
                           "dipg" = seuobj_dipg,
                           "mb" = seuobj_mb)

#identifying features
features <- SelectIntegrationFeatures(object.list = tumors.seuobj.list)
#creating anchors
tumors.anchors <- FindIntegrationAnchors(object.list = tumors.seuobj.list, anchor.features = features)
#creating integrated obj
tumors.seuobj <- IntegrateData(anchorset = tumors.anchors)


#running standard analysis workflow
tumors.seuobj <- NormalizeData(tumors.seuobj)
tumors.seuobj <- FindVariableFeatures(tumors.seuobj)
tumors.seuobj <- RunFastMNN(object.list = SplitObject(tumors.seuobj, split.by = "orig.ident"))
tumors.seuobj <- ScaleData(tumors.seuobj)
tumors.seuobj <- RunPCA(tumors.seuobj)
ElbowPlot(tumors.seuobj)
DimHeatmap(tumors.seuobj, dims = 1:15, cells = 500, balanced = TRUE)

tumors.seuobj <- RunUMAP(tumors.seuobj, dims = 1:15)
tumors.seuobj <- FindNeighbors(tumors.seuobj, dims = 1:15)
tumors.seuobj <- FindClusters(tumors.seuobj)
tumors.seuobj <- FindClusters(tumors.seuobj, resolution = 0.5)

DimPlot(tumors.seuobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(tumors.seuobj, reduction = "umap", group.by = "orig.ident")
DimPlot(tumors.seuobj, reduction = "umap", group.by = "mnn.reconstructed_snn_res.0.5")

#reset orig.ident
tumors.seuobj@meta.data$orig.ident <- gsub("SeuratProject", "MB", tumors.seuobj@meta.data$orig.ident)

#Differentially Expressed Genes in tumors.seuobj=============
tumors.seuobj.markers <- FindAllMarkers(tumors.seuobj)

tumors.seuobj.markers %>%
  group_by(cluster) %>% #grouping dataframe by cluster column
  dplyr::filter(avg_log2FC > 1) %>% #filtering all groups by log2FC > 1
  slice_head(n = 10) %>% #selecting top 10 genes for each cluster
  ungroup() -> top10.tumors #restacking each group (aka cluster) into one variable dataframe i.e. "top10"

DoHeatmap(tumors.seuobj, features = top10.tumors$gene) + NoLegend()

#Differentially Expressed Genes in day 115=======================
seuobj115.markers <- FindAllMarkers(seuobj115)

seuobj115.markers %>%
  group_by(cluster) %>% #grouping dataframe by cluster column
  dplyr::filter(avg_log2FC > 1) %>% #filtering all groups by log2FC > 1
  slice_head(n = 10) %>% #selecting top 10 genes for each cluster
  ungroup() -> top10.115 #restacking each group (aka cluster) into one variable dataframe i.e. "top10"

DoHeatmap(seuobj115, features = top10.115$gene) + NoLegend()

#gene profile for each cell type cluster in day 115
oligodendrocyte.genes <- top10.115$gene[1:10]
SLC24A4_PEX5L.genes <- top10.115$gene[11:20]
vascular_endothelial.genes <- top10.115$gene[21:30]
astrocyte.genes <- top10.115$gene[31:40]
granule_neuron.genes <- top10.115$gene[41:50]
purkinje_neuron.genes <- top10.115$gene[51:60]
inhibitory_interneuron.genes <- top10.115$gene[61:70]
unipolar_brush.genes <- top10.115$gene[71:80]
microglia.genes <- top10.115$gene[81:90]

#feature plot for each cell type gene profile on tumor cells
FeaturePlot(tumors.seuobj, features = oligodendrocyte.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = SLC24A4_PEX5L.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = vascular_endothelial.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = astrocyte.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = granule_neuron.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = purkinje_neuron.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = inhibitory_interneuron.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = unipolar_brush.genes, reduction = "umap", min.cutoff="q1")
FeaturePlot(tumors.seuobj, features = microglia.genes, reduction = "umap", min.cutoff="q1")


#gene profile for each cluster in tumor
cluster0.genes <- top10.tumors$gene[1:10]
cluster1.genes <- top10.tumors$gene[11:20]
cluster2.genes <- top10.tumors$gene[21:30]
cluster3.genes <- top10.tumors$gene[31:40]
cluster4.genes <- top10.tumors$gene[41:50]
cluster5.genes <- top10.tumors$gene[51:60]
cluster6.genes <- top10.tumors$gene[61:70]
cluster7.genes <- top10.tumors$gene[71:80]
cluster8.genes <- top10.tumors$gene[81:90]
cluster9.genes <- top10.tumors$gene[91:100]
cluster10.genes <- top10.tumors$gene[101:110]
cluster11.genes <- top10.tumors$gene[111:120]


#Adding Module Scores for each cell type onto tumor cells=============
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=purkinje_neuron.genes), assay = "RNA", name = "purkinje_neuron")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=oligodendrocyte.genes), assay = "RNA", name = "oligodendrocyte")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=SLC24A4_PEX5L.genes), assay = "RNA", name = "SLC24A4_PEX5L")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=vascular_endothelial.genes), assay = "RNA", name = "vascular_endothelial")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=astrocyte.genes), assay = "RNA", name = "astrocyte")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=granule_neuron.genes), assay = "RNA", name = "granule_neuron")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=inhibitory_interneuron.genes), assay = "RNA", name = "inhibitory_interneuron")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=unipolar_brush.genes), assay = "RNA", name = "unipolar_brush")
tumors.seuobj <- AddModuleScore(tumors.seuobj, features = list("celltype"=microglia.genes), assay = "RNA", name = "microglia")


##visualization of module scores
#Dot Plot
DotPlot(tumors.seuobj, features = c("purkinje_neuron", "oligodendrocyte", "SLC24A4_PEX5L", "vascular_endothelial", "astrocyte", "granule_neuron", "inhibitory_interneuron", "unipolar_brush", "microglia"), group.by = "mnn.reconstructed_snn_res.0.8") + RotatedAxis()

cell_types <- c("purkinje_neuron", "oligodendrocyte", "SLC24A4_PEX5L", "vascular_endothelial", "astrocyte", "granule_neuron", "inhibitory_interneuron", "unipolar_brush", "microglia")

#Violin Plot
for(cell_type in cell_types){
  print(VlnPlot(tumors.seuobj, features = c(cell_type), group.by = "mnn.reconstructed_snn_res.0.8"))
}

#Feature Plot
for(cell_type in cell_types){
  print(FeaturePlot(tumors.seuobj, features = c(cell_type), min.cutoff = "q1"))
}


#Subset by cancer type======================

tumors.seuobj.mb <- subset(tumors.seuobj, subset = orig.ident == "MB")
tumors.seuobj.atrt <- subset(tumors.seuobj, subset = orig.ident == "atrt40" | orig.ident == "atrt41" | orig.ident == "atrt50" | orig.ident == "atrt61")
tumors.seuobj.dipg <- subset(tumors.seuobj, subset = orig.ident == "Glioma")

#Violin plots for medulloblastoma
for(cell_type in cell_types){
  print(VlnPlot(tumors.seuobj.mb, features = c(cell_type), group.by = "mnn.reconstructed_snn_res.0.8"))
}

#Violin plots for ATRT
for(cell_type in cell_types){
  print(VlnPlot(tumors.seuobj.atrt, features = c(cell_type), group.by = "mnn.reconstructed_snn_res.0.8"))
}

#Violin plots for DIPG
for(cell_type in cell_types){
  print(VlnPlot(tumors.seuobj.dipg, features = c(cell_type), group.by = "mnn.reconstructed_snn_res.0.8"))
}


#DimPlots by cluster
DimPlot(tumors.seuobj.mb, group.by = "mnn.reconstructed_snn_res.0.8", reduction = "umap")
DimPlot(tumors.seuobj.atrt, group.by = "mnn.reconstructed_snn_res.0.8", reduction = "umap")
DimPlot(tumors.seuobj.dipg, group.by = "mnn.reconstructed_snn_res.0.8", reduction = "umap")