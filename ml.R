# load libraries-----------------
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(ggplot2)

install.packages("ggpmisc")

library(ggpmisc)
library(nlme)

basepath <- "/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/LogOdds_Plots_Astrocytes/"

# create integrated object

seuobj115 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj115.RDS")
seuobj125 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj125.RDS")

subsets.seuobj.list <- list("115" = seuobj115,
                            "125" = seuobj125)

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

#creating data frame of gene expression data
day115_125_gene_expression_data <- as.data.frame(subsets.seuobj@assays$integrated$data)
# switching cols and rows by transposing
day115_125_gene_expression_data <- t(day115_125_gene_expression_data)
# make data remain as data.frame
day115_125_gene_expression_data <- data.frame(day115_125_gene_expression_data)
# round vals to 2 spots
day115_125_gene_expression_data <- (round(day115_125_gene_expression_data, 2))
# add cell id col to metadata of seuobj
subsets.seuobj$cellid <- rownames(subsets.seuobj@meta.data)
# making a metadata dataframe of seuobj with main cluster and cell id
metadata_df <- subsets.seuobj@meta.data[, c("Main_cluster_name", "cellid")]
# adding cell ids to metadata
day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)
# combining data and metadata
day115_125_gene_expression_data <- cbind(day115_125_gene_expression_data, metadata_df)
# removing the cell id from the expression dataframe
day115_125_gene_expression_data$cellid <- NULL

# log plots
count <- unique(day115_125_gene_expression_data[,1])
count_df <- data.frame(count=count, percent=rep(NA, length(count)), logodds=rep(NA,length(count)))


for (gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"][1:101]) {
  # n_astro_list <- c()
  message(gene)
  count <- unique(day115_125_gene_expression_data[,gene])
  
  count_df = data.frame(count=count, 
                        percent = rep(NA, length(count)), 
                        logodds= rep(NA, length(count)))
  
  message("Count_df created")
  for(n in unique(count_df$count)){
    # count for all astrocytes
    n_astrocytes <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Astrocytes",gene]==n)
  
    # probability of a cell being an astrocyte
    # number of astrocytes w unique count value for spp1 / all cells
    percent <- n_astrocytes/ sum(day115_125_gene_expression_data[,gene]==n)
  
    # log odds
    logodds <- log(percent/(1-percent))
  
    # adding all variables to count_df
    count_df[count_df$count == n, "percent"] <- percent
    count_df[count_df$count == n, "logodds"] <- logodds
  }
  # print(count_df)
  # table(count_df$logodds==-Inf)

  count_df <- count_df[count_df$logodds != -Inf, ]
  
  if(nrow(count_df) >= 20) {
    plot <- ggplot(count_df, aes(x = count, y = logodds)) + 
      geom_point(stat = "identity") + 
      geom_smooth(method = "lm") + 
      ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) +
      stat_poly_eq(use_label(c("eq", "R2")))
    
  # saveRDS(plot, paste0(basepath, gene, ".RDS"))
    saveRDS(plot, paste0(basepath, "plot.RDS"))
    
    png(paste0(basepath, gene, ".png"), width=300, height=300)
    print(plot)
    dev.off()
  } 
}