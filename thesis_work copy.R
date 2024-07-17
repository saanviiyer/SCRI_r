# load libraries-----------------
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(sctransform)

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

seuobj_atrt@active.assay = "RNA"
seuobj_dipg@active.assay = "RNA"
seuobj_mb@active.assay = "RNA"
seuobj89@active.assay = "RNA"
seuobj94@active.assay = "RNA"
seuobj110@active.assay = "RNA"
seuobj115@active.assay = "RNA"
seuobj125@active.assay = "RNA"


#atrt pre processing
seuobj_atrt[["percent.mt"]] <- PercentageFeatureSet(seuobj_atrt, pattern = "^MT-")
# seuobj_atrt <- subset(seuobj_atrt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seuobj_atrt <- NormalizeData(seuobj_atrt, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_atrt <- FindVariableFeatures(seuobj_atrt, selection.method = "vst", nfeatures = 2000)
seuobj_atrt <- ScaleData(seuobj_atrt, features = rownames(seuobj_atrt))
seuobj_atrt <- RunPCA(seuobj_atrt, features = VariableFeatures(object = seuobj_atrt))
seuobj_atrt <- FindNeighbors(seuobj_atrt, dims = 1:15)
seuobj_atrt <- FindClusters(seuobj_atrt, resolution = 0.8)
seuobj_atrt <- RunUMAP(seuobj_atrt, dims = 1:15)

#dipg pre processing
seuobj_dipg[["percent.mt"]] <- PercentageFeatureSet(seuobj_dipg, pattern = "^MT-")
seuobj_dipg <- subset(seuobj_dipg, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seuobj_dipg <- NormalizeData(seuobj_dipg, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_dipg <- FindVariableFeatures(seuobj_dipg, selection.method = "vst", nfeatures = 2000)
seuobj_dipg <- ScaleData(seuobj_dipg, features = rownames(seuobj_dipg))
seuobj_dipg <- RunPCA(seuobj_dipg, features = VariableFeatures(object = seuobj_dipg))
seuobj_dipg <- FindNeighbors(seuobj_dipg, dims = 1:15)
seuobj_dipg <- FindClusters(seuobj_dipg, resolution = 0.8)
seuobj_dipg <- RunUMAP(seuobj_dipg, dims = 1:15)

#mb pre processing
seuobj_mb <- NormalizeData(seuobj_mb, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_mb <- FindVariableFeatures(seuobj_mb, selection.method = "vst", nfeatures = 2000)
seuobj_mb <- ScaleData(seuobj_mb, features = rownames(seuobj_mb))
seuobj_mb <- RunPCA(seuobj_mb, features = VariableFeatures(object = seuobj_mb))
seuobj_mb <- FindNeighbors(seuobj_mb, dims = 1:15)
seuobj_mb <- FindClusters(seuobj_mb, resolution = 0.8)
seuobj_mb <- RunUMAP(seuobj_mb, dims = 1:15)

# sc transform
seuobj_atrt <- SCTransform(seuobj_atrt, vars.to.regress = "percent.mt", verbose = FALSE)
seuobj_dipg <- SCTransform(seuobj_dipg, vars.to.regress = "percent.mt", verbose = FALSE)

#Integration of tumors===========
#creating tumor seuobj list
seuobj.list <- list("atrt" = seuobj_atrt,
                           "dipg" = seuobj_dipg,
                           "mb" = seuobj_mb,
                           "89" = seuobj89,
                           "94" = seuobj94,
                           "110" = seuobj110,
                           "115" = seuobj115,
                           "125" = seuobj125
                           )

seuobj_atrt@active.assay = "SCT"
seuobj_dipg@active.assay = "SCT"
seuobj_mb@active.assay = "SCT"
seuobj89@active.assay = "SCT"
seuobj94@active.assay = "SCT"
seuobj110@active.assay = "SCT"
seuobj115@active.assay = "SCT"
seuobj125@active.assay = "SCT"

#identifying features
features <- SelectIntegrationFeatures(object.list = seuobj.list)
#creating anchors
tumors.anchors <- FindIntegrationAnchors(object.list = seuobj.list, anchor.features = features)
#creating integrated obj
tumors.seuobj <- IntegrateData(anchorset = tumors.anchors)

tumors.seuobj <- NormalizeData(tumors.seuobj, normalization.method = "LogNormalize", scale.factor = 10000)
tumors.seuobj <- FindVariableFeatures(tumors.seuobj, selection.method = "vst", nfeatures = 2000)
tumors.seuobj <- RunFastMNN(object.list = SplitObject(tumors.seuobj, split.by = "orig.ident"))

tumors.seuobj <- ScaleData(tumors.seuobj, features = rownames(tumors.seuobj))
tumors.seuobj <- RunPCA(tumors.seuobj, features = VariableFeatures(object = tumors.seuobj))
tumors.seuobj <- FindNeighbors(tumors.seuobj, dims = 1:15)
tumors.seuobj <- FindClusters(tumors.seuobj, resolution = 0.8)
tumors.seuobj <- RunUMAP(tumors.seuobj, dims = 1:15)

# DimPlot(tumors.seuobj, group.by = "Main_cluster_name", label = TRUE)
DimPlot(tumors.seuobj, group.by = "orig.ident", label = TRUE)

integrated.seuobj <- tumors.seuobj

basepath <- "/Users/saanviiyer/Downloads/org_data/"
saveRDS(tumors.seuobj, paste0(basepath, "integrated_seuratobj.RDS"))

# seuobj125@meta.data$Topic__33

# FeaturePlots
FeaturePlot(seuobj125, features = "Topic__33", reduction = "umap")
FeaturePlot(seuobj115, features = "Topic__12", reduction = "umap")
FeaturePlot(seuobj110, features = "Topic__14", reduction = "umap")
FeaturePlot(seuobj94, features = "Topic__27", reduction = "umap")
FeaturePlot(seuobj89, features = "Topic__14", reduction = "umap")

table(integrated.seuobj@meta.data$seurat_clusters)

DimPlot(integrated.seuobj, group.by="seurat_clusters", label=TRUE)

# DEG
integrated.seuobj <- FindAllMarkers(integrated.seuobj)
i.s <- integrated.seuobj
i.s_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
