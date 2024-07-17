# load libraries-----------------
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)

# load data----------------------
#developmental days
seuobj94 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj94.RDS")
seuobj94_astrocyte <- subset(seuobj94, subset = Main_cluster_name == 'Astrocytes')
seuobj94_granule <- subset(seuobj94, subset = Main_cluster_name == 'Granule neurons')

seuobj_atrt <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_atrt.RDS")
seuobj_dipg <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_dipg.RDS")
seuobj_mb <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_mb.RDS")

subsets.seuobj.list <- list("ast" = seuobj94_astrocyte,
                           "gra" = seuobj94_granule)

#identifying features
features <- SelectIntegrationFeatures(object.list = subsets.seuobj.list)
#creating anchors
subsets.anchors <- FindIntegrationAnchors(object.list = subsets.seuobj.list, anchor.features = features)
#creating integrated obj
subsets.seuobj <- IntegrateData(anchorset = subsets.anchors)

subsets.seuobj <- NormalizeData(subsets.seuobj)
subsets.seuobj <- FindVariableFeatures(subsets.seuobj, selection.method = "vst", nfeatures = 2000)
#Scaling---------------------------------------------------
all.genes <- rownames(subsets.seuobj)
subsets.seuobj <- ScaleData(subsets.seuobj, features = all.genes)
#PCA--------------------------------------------------------
subsets.seuobj <- RunPCA(subsets.seuobj, features = VariableFeatures(object = subsets.seuobj))

PC1_positive = c("CST3", "KCND2", "ERBB4", "VCAM1", "CADM2", "AQP4", "NKAIN3", "SLC6A11", "GJA1", "NTRK2",
"NRXN1", "NCAN", "BCAN", "HSPG2", "KCND3", "TFPI", "RNF219-AS1", "PPAP2B", "GLIS3", "PI15",
"CD44", "SPON1", "INHBB", "SGCD", "RP11-689K5.3", "GPC5", "ATP1B2", "DOCK10", "CPXM2", "LRP1B")

PC1_negative = c("LRRTM4", "DACH1", "RP11-649A16.1", "DGKB", "GRID2", "NRG1", "CTD-3088G3.8", "TNC", "ESRRG", "PRTG",
                 "SHROOM3", "CDH4", "PTPRM", "HS3ST4", "NELL2", "GRM3", "PCSK6", "PRKG1", "CTNNA3", "CDH10",
                 "COL23A1", "FREM2", "RERG", "NCKAP5", "NRP1", "THSD4", "SLC1A3", "ADAMTS9", "RALYL", "DMD")
PC2_positive = c("NRXN1", "KCND2", "ERBB4", "ANKFN1", "GRIK1", "RNF219-AS1", "ANKS1B", "EGFR", "DCLK2", "MEGF11",
                 "SNCAIP", "RFTN2", "LUZP2", "NKAIN3", "CSMD1", "TMTC2", "LSAMP", "ZNF521", "ZBTB20", "KCNH7",
                 "SOX2-OT", "LINC00478", "GLIS3", "RYR3", "ADAMTSL1", "RP11-553L6.5", "DPF3", "LRRTM3", "LRP1B", "NHSL1")
PC2_negative = c("HSPG2", "CD44", "TMBIM1", "VCAM1", "MEPCE", "RP11-642D21.2", "CETN2", "CPLX4", "NFATC4", "DDIT4",
                 "S100B", "HSPA1B", "TMEM205", "NCDN", "HIST1H2BK", "FAM107A", "METTL7A", "IL6R", "S1PR3", "FYCO1",
                 "PLCD1", "TSTA3", "DNPH1", "RP11-815J21.4", "LECT1", "ALYREF", "NOP10", "OLFM3", "ARHGDIA", "PEF1")
PC5_positive = c("COL4A2", "TNFRSF1A", "KCNJ16", "TMEM115", "CNKSR3", "SFRP2", "NXPH1", "TMED1", "TNFRSF11B", "TTC29",
                 "RP11-630C16.2", "MERTK", "ADAMTS17", "FGFRL1", "FIBIN", "RNU7-40P", "SLC7A11", "MKI67IP", "WWTR1", "RILPL2",
                 "RP11-572C15.6", "SLC35A2", "RP11-522B15.3", "RP11-133F8.2", "CPVL", "TM2D2", "ANGPTL2", "TTC6", "TMEM144", "C10orf11")

PC5_negative = c("NACC2", "MAP3K19", "MASP1", "DNAH7", "KCNJ10", "LAMB2", "KIAA1161", "RPS6", "NNAT", "GRIA3",
                 "RAVER1", "ATP10B", "TP53I13", "SEC61G", "INPP5E", "FAM198B", "ZFP36L1", "RP11-436K8.1", "TMEM205", "EEPD1",
                 "HYAL2", "AC002429.5", "UST", "VEGFC", "NDNL2", "KCTD12", "SEMA5B", "ECM2", "MEGF10", "INTS5")


subsets.seuobj <- AddModuleScore(subsets.seuobj, features = list("PC"=PC1_positive), assay = "RNA", name = "PC1_positive")
subsets.seuobj <- AddModuleScore(subsets.seuobj, features = list("PC"=PC1_negative), assay = "RNA", name = "PC1_negative")

subsets.seuobj <- AddModuleScore(subsets.seuobj, features = list("PC"=PC2_positive), assay = "RNA", name = "PC2_positive")
subsets.seuobj <- AddModuleScore(subsets.seuobj, features = list("PC"=PC2_negative), assay = "RNA", name = "PC2_negative")

subsets.seuobj <- AddModuleScore(subsets.seuobj, features = list("PC"=PC2_positive), assay = "RNA", name = "PC5_positive")
subsets.seuobj <- AddModuleScore(subsets.seuobj, features = list("PC"=PC2_negative), assay = "RNA", name = "PC5_negative")

pcs <- c("PC1_positive1", "PC1_negative1", "PC2_positive1", "PC2_negative1", "PC5_positive1", "PC5_negative1")
subsets.seuobj <- RunUMAP(subsets.seuobj, dims = 1:10)  # Adjust dimensions as necessary

#Feature Plot
for(pc in pcs){
  print(FeaturePlot(subsets.seuobj, features = c(pc), reduction = "umap", min.cutoff = "q1"))
}


DimPlot(subsets.seuobj, group.by = "Main_cluster_name", label = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------

# load libraries-----------------
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)

# load data----------------------

#tumors
seuobj_atrt <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_atrt.RDS")
seuobj_dipg <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_dipg.RDS")
seuobj_mb <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj_mb.RDS")

#atrt pre processing
seuobj_atrt[["percent.mt"]] <- PercentageFeatureSet(seuobj_atrt, pattern = "^MT-")
# seuobj_atrt <- subset(seuobj_atrt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
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
seuobj_dipg@active.assay="RNA"
seuobj_dipg <- NormalizeData(seuobj_dipg, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_dipg <- FindVariableFeatures(seuobj_dipg, selection.method = "vst", nfeatures = 2000)
seuobj_dipg <- ScaleData(seuobj_dipg, features = rownames(seuobj_dipg))
seuobj_dipg <- RunPCA(seuobj_dipg, features = VariableFeatures(object = seuobj_dipg))

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

tumors.seuobj <- RunUMAP(tumors.seuobj, dims = 1:15)
tumors.seuobj <- FindNeighbors(tumors.seuobj, dims = 1:15)
tumors.seuobj <- FindClusters(tumors.seuobj)
tumors.seuobj <- FindClusters(tumors.seuobj, resolution = 0.5)

DimPlot(tumors.seuobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(tumors.seuobj, reduction = "umap", group.by = "orig.ident")
DimPlot(tumors.seuobj, reduction = "umap", group.by = "mnn.reconstructed_snn_res.0.5")


