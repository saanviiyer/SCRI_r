library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(sctransform)

# load data----------------------
#developmental days
seuobj_cluster<- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/seuobj0_17_18_20_23.RDS")

Day89_Topic_14_top15 <- c("NLGN1", "LSAMP", "ADAMTS6", "NFIB", "NFIA", "LINC00669", "RP11-649A16.1", "SYNE2", "EGFEM1P", "DIAPH3", "DLEU2", "WWOX", "CENPF", "ZBTB20", "RAD51B")
Day94_Topic_27_top15 <- c("SYNE2", "LINC00669", "EGFEM1P", "DIAPH3", "CENPF", "DLEU2", "MKI67", "FBXL7", "NFIA", "CENPP", "NLGN1", "TOP2A", "CEP128", "RFC3", "ASPM")
Day110_Topic14_top15 <- c("LINC00669", "SYNE2", "EGFEM1P", "NLGN1", "DIAPH3", "CENPP", "DLEU2", "WWOX", "RFC3", "APOLD1", "TOP2A", "CENPF", "BRIP1", "FBXL7", "MKI67")
Day115_Topic_12 <- c("LRRTM4", "PTN", "DACH1", "SLC1A3", "TNC", "PTPRM", "GPM6B", "RP11-649A16.1", "SEMA6A", "CTNNA3", "EGFEM1P", "CDH4", "NCKAP5", "NAALADL2", "DGKB")
Day125_Topic_33 <- c("LINC00669", "LSAMP", "NLGN1", "DCC", "DIAPH3", "CENPP", "DLEU2", "RFC3", "MKI67", "APOLD1", "WWOX", "TOP2A", "CENPF", "RAD51B", "EGFEM1P")

top15s <- c(Day89_Topic_14_top15, Day94_Topic_27_top15, Day110_Topic14_top15, Day115_Topic_12, Day125_Topic_33)
top15names <- c("Day89_Topic_14_top151", "Day94_Topic_27_top151","Day110_Topic14_top151", "Day115_Topic_121", "Day125_Topic_331")

for(i in 1:5) {
  seuobj_cluster <- AddModuleScore(seuobj_cluster, features = list("genes" = top15s[i]), assay = "RNA", name = top15names[i])
}

for(j in top15names){
  print(FeaturePlot(seuobj_cluster, features = c(j), reduction="umap"))
  print(  DimPlot(seuobj_cluster,group.by="Main_cluster_name"))
  print(  DimPlot(seuobj_cluster,group.by="sub_cluster_name"))
}

DimPlot(seuobj_cluster,group.by="orig.ident")