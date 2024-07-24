# load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(ggplot2)

install.packages("ggpmisc")
install.packages("broom")
library(broom)
library(ggpmisc)
library(nlme)

basepath <- "/Users/saanviiyer/Documents/GitHub/SCRI_r/LogOdds_Plots_Astrocytes/"

# create integrated object

seuobj115 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj115.RDS")
seuobj125 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj125.RDS")
seuobj110 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj110.RDS")

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
#Scaling-
all.genes <- rownames(subsets.seuobj)
subsets.seuobj <- ScaleData(subsets.seuobj, features = all.genes)
#PCA
subsets.seuobj <- RunPCA(subsets.seuobj, features = VariableFeatures(object = subsets.seuobj))


seuobj110 <- NormalizeData(seuobj110)
all.genes <- rownames(seuobj110)
seuobj110 <- ScaleData(seuobj110, features = all.genes)
#PCA
seuobj110 <- RunPCA(seuobj110, features = VariableFeatures(object = subsets.seuobj))
day110_gene_expression_data <- as.data.frame(seuobj110@assays$RNA$data)


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

rsq_df <- data.frame(gene=colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"],
                     rsquare=rep(NA, 2000))

for(gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"]){
  
  message(gene)
  
  count <- unique(day115_125_gene_expression_data[,gene])
  
  count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
  
  message('Count df created')
  
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
  count_df <- count_df[count_df$logodds != Inf, ]
  if(nrow(count_df) >= 20){
    
    
    model <- lm(logodds ~ count, data = count_df)
    model_summary <- glance(model)
    r_squared <- model_summary$r.squared
    
    rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared
    
    # plot <- ggplot(count_df, aes(x = count, y = logodds)) +
    #   geom_point(stat = "identity") + geom_smooth(method = "lm") +
    #   ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
    # 
    # 
    # 
    # png(paste0(basepath, gene, ".png"), width = 300, height = 300)
    # print(plot)
    # dev.off()
  }
}


rsq_df <- rsq_df[order(-rsq_df$rsquare),]

# for(gene in rsq_df$gene[1:30]){
for(i in 1 : 30) {

  if (rsq_df$rsquare[i] >= 0.4) {
    print(paste0("COUNT",i))
    gene <- rsq_df$gene[i]
    message(gene)
    
    count <- unique(day115_125_gene_expression_data[,gene])
    
    count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
    
    message('Count df created')
    
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
    count_df <- count_df[count_df$logodds != Inf, ]
    
    if(nrow(count_df) >= 20){
      
      
      # model <- lm(logodds ~ count, data = count_df)
      # model_summary <- glance(model)
      # r_squared <- model_summary$r.squared
      # 
      # rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared
      
      plot <- ggplot(count_df, aes(x = count, y = logodds)) +
        geom_point(stat = "identity") + geom_smooth(method = "lm") +
        ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
      
      
      
      png(paste0(basepath, gene, "TOP30",".png"), width = 300, height = 300)
      print(plot)
      dev.off()
    }
  }
}
# 
# # creating top_genes_logodds
# top_genes_logodds <- c()
# 
# # populating top_genes_logodds with top genes
# for(i in 1 : 30) {
# 
#   if (rsq_df$rsquare[i] >= 0.2) {
#     top_genes_logodds <- append(top_genes_logodds, rsq_df$gene[i])
#   }
# }
# 
# day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)
# 
# # split 75% training data
# training.data <- day115_125_gene_expression_data %>%
#   group_by(Main_cluster_name) %>%
#   filter(row_number() <= 0.75 * n())
# 
# # all data not in training wil lbe used for testing
# testing.data <- day115_125_gene_expression_data[!day115_125_gene_expression_data$cellid %in% training.data$cellid,]
# 
# rownames(training.data) <- training.data$cellid
# rownames(testing.data) <- testing.data$cellid
# 
# training.data$cellid <- NULL
# testing.data$cellid <- NULL
# 
# # creating top genes from each
# training.data.topgenes <- training.data[,c(top_genes_logodds, "Main_cluster_name")]
# testing.data.topgenes <- testing.data[,c(top_genes_logodds, "Main_cluster_name")]
# 
# # training.data.topgenes$Main_cluster_name <- 
# 
# training.data.topgenes$celltype <- 0
# training.data.topgenes$celltype[training.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1
# 
# testing.data.topgenes$celltype <- 0
# testing.data.topgenes$celltype[testing.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1
# 
# training.data.topgenes$Main_cluster_name = NULL
# testing.data.topgenes$Main_cluster_name = NULL
# 
# # training GLM with 0.75 GE data
# model_lm <- glm(celltype ~. , data = training.data.topgenes, family = "binomial")
# 
# # testing GLM with 0.25 GE data
# predictions <- model_lm %>% predict(testing.data.topgenes)
# 
# # determining the model performance
# performance <- data.frame(
#   RMSE = RMSE(predictions, testing.data.topgenes$celltype),
#   R2 = R2(predictions, testing.data.topgenes$celltype)
# )
#   
# summary(model_lm)  

#getting the top genes from rsq_)df
topgenes_logodds <- c()
for(i in 1 : 30) {
  if (rsq_df$rsquare[i] >= 0.4) {
    topgenes_logodds <- append(topgenes_logodds, rsq_df$gene[i])
  }
}

newgenes <- c("S100B", "GFAP", "APOE", "AQP4")

topgenes_logodds <- append(topgenes_logodds, newgenes)

#load dplyr
library(dplyr)

#making the cell id columns the same as the rownames
day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)

#making the training dat groping by main clusrter name and filtering out
training.data <- day115_125_gene_expression_data %>%
  group_by(Main_cluster_name) %>%
  filter(row_number() <= 0.75 * n())

#making testing data from everything that isn't training
testing.data <- day115_125_gene_expression_data[!day115_125_gene_expression_data$cellid %in% training.data$cellid, ]

rownames(training.data) <- training.data$cellid
rownames(testing.data) <- testing.data$cellid

training.data$cellid <- NULL
testing.data$cellid <- NULL


#creating top genes from each
training.data.topgenes <- training.data[, c(topgenes_logodds, "Main_cluster_name")]
testing.data.topgenes <- testing.data[, c(topgenes_logodds, "Main_cluster_name")]

#creating binary predictor column (astrocyte = 1, other = 0)
training.data.topgenes$celltype <- 0
training.data.topgenes$celltype[training.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1

testing.data.topgenes$celltype <- 0
testing.data.topgenes$celltype[testing.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1

training.data.topgenes$Main_cluster_name <- NULL
testing.data.topgenes$Main_cluster_name <- NULL


#training glm with 0.75 GE data
model_lm <- glm(celltype ~. , data = training.data.topgenes, family = "binomial")

#testing GLM with 0.25 GE data
predictions <- model_lm %>% predict(testing.data.topgenes)

#determining the model performance
performance <- data.frame(RMSE = RMSE(predictions, testing.data.topgenes$celltype),
                          R2 = R2(predictions, testing.data.topgenes$celltype))

#view model summary and performance
summary(model_lm)
performance

# VIF to determine multicolinearity with car package
library(car)

car::vif(model_lm)

# stepwise AIC
library(MASS)
# stepwise regression --forward
step.model.forward <- stepAIC(model_lm, direction = "forward", trace = 10)
summary(step.model.forward)

# backwards
step.model.backward <- stepAIC(model_lm, direction = "backward", trace = 10)
summary(step.model.backward)

# both
step.model.both <- stepAIC(model_lm, direction = "both", trace = 10)
summary(step.model.both)

# NEURONAL MODEL ---------------------


# log plots
count <- unique(day115_125_gene_expression_data[,1])
count_df <- data.frame(count=count, percent=rep(NA, length(count)), logodds=rep(NA,length(count)))

rsq_df <- data.frame(gene=colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"],
                     rsquare=rep(NA, 2000))

for(gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"]){
  
  message(gene)
  
  count <- unique(day115_125_gene_expression_data[,gene])
  
  count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
  
  message('Count df created')
  
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
  count_df <- count_df[count_df$logodds != Inf, ]
  if(nrow(count_df) >= 20){
    
    
    model <- lm(logodds ~ count, data = count_df)
    model_summary <- glance(model)
    r_squared <- model_summary$r.squared
    
    rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared
    
    # plot <- ggplot(count_df, aes(x = count, y = logodds)) +
    #   geom_point(stat = "identity") + geom_smooth(method = "lm") +
    #   ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
    # 
    # 
    # 
    # png(paste0(basepath, gene, ".png"), width = 300, height = 300)
    # print(plot)
    # dev.off()
  }
}


rsq_df <- rsq_df[order(-rsq_df$rsquare),]

# for(gene in rsq_df$gene[1:30]){
for(i in 1 : 30) {
  
  if (rsq_df$rsquare[i] >= 0.4) {
    print(paste0("COUNT",i))
    gene <- rsq_df$gene[i]
    message(gene)
    
    count <- unique(day115_125_gene_expression_data[,gene])
    
    count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
    
    message('Count df created')
    
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
    count_df <- count_df[count_df$logodds != Inf, ]
    
    if(nrow(count_df) >= 20){
      
      
      # model <- lm(logodds ~ count, data = count_df)
      # model_summary <- glance(model)
      # r_squared <- model_summary$r.squared
      # 
      # rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared
      
      plot <- ggplot(count_df, aes(x = count, y = logodds)) +
        geom_point(stat = "identity") + geom_smooth(method = "lm") +
        ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
      
      
      
      png(paste0(basepath, gene, "TOP30",".png"), width = 300, height = 300)
      print(plot)
      dev.off()
    }
  }
}
# 
# # creating top_genes_logodds
# top_genes_logodds <- c()
# 
# # populating top_genes_logodds with top genes
# for(i in 1 : 30) {
# 
#   if (rsq_df$rsquare[i] >= 0.2) {
#     top_genes_logodds <- append(top_genes_logodds, rsq_df$gene[i])
#   }
# }
# 
# day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)
# 
# # split 75% training data
# training.data <- day115_125_gene_expression_data %>%
#   group_by(Main_cluster_name) %>%
#   filter(row_number() <= 0.75 * n())
# 
# # all data not in training wil lbe used for testing
# testing.data <- day115_125_gene_expression_data[!day115_125_gene_expression_data$cellid %in% training.data$cellid,]
# 
# rownames(training.data) <- training.data$cellid
# rownames(testing.data) <- testing.data$cellid
# 
# training.data$cellid <- NULL
# testing.data$cellid <- NULL
# 
# # creating top genes from each
# training.data.topgenes <- training.data[,c(top_genes_logodds, "Main_cluster_name")]
# testing.data.topgenes <- testing.data[,c(top_genes_logodds, "Main_cluster_name")]
# 
# # training.data.topgenes$Main_cluster_name <- 
# 
# training.data.topgenes$celltype <- 0
# training.data.topgenes$celltype[training.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1
# 
# testing.data.topgenes$celltype <- 0
# testing.data.topgenes$celltype[testing.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1
# 
# training.data.topgenes$Main_cluster_name = NULL
# testing.data.topgenes$Main_cluster_name = NULL
# 
# # training GLM with 0.75 GE data
# model_lm <- glm(celltype ~. , data = training.data.topgenes, family = "binomial")
# 
# # testing GLM with 0.25 GE data
# predictions <- model_lm %>% predict(testing.data.topgenes)
# 
# # determining the model performance
# performance <- data.frame(
#   RMSE = RMSE(predictions, testing.data.topgenes$celltype),
#   R2 = R2(predictions, testing.data.topgenes$celltype)
# )
#   
# summary(model_lm)  

#getting the top genes from rsq_)df
topgenes_logodds <- c()
for(i in 1 : 30) {
  if (rsq_df$rsquare[i] >= 0.4) {
    topgenes_logodds <- append(topgenes_logodds, rsq_df$gene[i])
  }
}

newgenes <- c("S100B", "GFAP", "APOE", "AQP4")

topgenes_logodds <- append(topgenes_logodds, newgenes)

#load dplyr
library(dplyr)

#making the cell id columns the same as the rownames
day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)

#making the training dat groping by main clusrter name and filtering out
training.data <- day115_125_gene_expression_data %>%
  group_by(Main_cluster_name) %>%
  filter(row_number() <= 0.75 * n())

#making testing data from everything that isn't training
testing.data <- day115_125_gene_expression_data[!day115_125_gene_expression_data$cellid %in% training.data$cellid, ]

rownames(training.data) <- training.data$cellid
rownames(testing.data) <- testing.data$cellid

training.data$cellid <- NULL
testing.data$cellid <- NULL


#creating top genes from each
training.data.topgenes <- training.data[, c(topgenes_logodds, "Main_cluster_name")]
testing.data.topgenes <- testing.data[, c(topgenes_logodds, "Main_cluster_name")]

#creating binary predictor column (astrocyte = 1, other = 0)
training.data.topgenes$celltype <- 0
training.data.topgenes$celltype[training.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1

testing.data.topgenes$celltype <- 0
testing.data.topgenes$celltype[testing.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1

training.data.topgenes$Main_cluster_name <- NULL
testing.data.topgenes$Main_cluster_name <- NULL


#training glm with 0.75 GE data
model_lm <- glm(celltype ~. , data = training.data.topgenes, family = "binomial")

#testing GLM with 0.25 GE data
predictions <- model_lm %>% predict(testing.data.topgenes)

#determining the model performance
performance <- data.frame(RMSE = RMSE(predictions, testing.data.topgenes$celltype),
                          R2 = R2(predictions, testing.data.topgenes$celltype))

#view model summary and performance
summary(model_lm)
performance

# VIF to determine multicolinearity with car package
library(car)

car::vif(model_lm)

# stepwise AIC
library(MASS)
# stepwise regression --forward
step.model.forward <- stepAIC(model_lm, direction = "forward", trace = 10)
summary(step.model.forward)

# backwards
step.model.backward <- stepAIC(model_lm, direction = "backward", trace = 10)
summary(step.model.backward)

# both
step.model.both <- stepAIC(model_lm, direction = "both", trace = 10)
summary(step.model.both)

----------------------------------
#   
#   #workflow for neurons
#   
# rsq_df_neurons <- data.frame(gene = colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"], rsquare = rep(NA, 2000))
# 
# rsq_df_neurons[]
# 
# 
# 
# 
# 
# for(gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"]){
#   #tells name of current gene
#   message(gene)
#   #filling count with gene expression values for all cells for this gene
#   count <- unique(day115_125_gene_expression_data[,gene])
#   #creating a dataframe with percent and logodds columns corresponding to probability that cell is a granule neuron or inhibitory interneuron depending on counts for that gene
#   count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
#   
#   #displays that coutn df was created
#   message('Count df created')
#   
#   #for loop iterating through each unique count value
#   for(n in unique(count_df$count)){
#     # count for all neurons
#     n_ii <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Inhibitory interneurons",gene]==n)
#     n_gn <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Granule neurons",gene]==n)
#     
#     n_neurons <- n_ii + n_gn
#     
#     # probability of a cell being an neuron
#     # number of neurons w unique count value for spp1 / all cells
#     percent <- n_neurons/ sum(day115_125_gene_expression_data[,gene]==n)
#     
#     # log odds
#     logodds <- log(percent/(1-percent))
#     
#     # adding all variables to count_df
#     count_df[count_df$count == n, "percent"] <- percent
#     count_df[count_df$count == n, "logodds"] <- logodds
#   }
#   # print(count_df)
#   # table(count_df$logodds==-Inf)
#   
#   #removing all -Inf logodds
#   count_df <- count_df[count_df$logodds != -Inf, ]
#   #removing all Inf logodds
#   count_df <- count_df[count_df$logodds != Inf, ]
#   
#   #only executing if there are more than 20 unique counts for the gene
#   if(nrow(count_df) >= 20){
#     
#     #creating a linear model to model relation between count and logodds
#     model <- lm(logodds ~ count, data = count_df)
#     #getting summary of the model
#     model_summary <- glance(model)
#     #extracting r squared from the summary
#     r_squared <- model_summary$r.squared
#     
#     #adding r62 value from the model into rsq_df
#     rsq_df_neurons[rsq_df_neurons$gene==gene, "rsquare"] <- r_squared
#     
#     # plot <- ggplot(count_df, aes(x = count, y = logodds)) +
#     #   geom_point(stat = "identity") + geom_smooth(method = "lm") +
#     #   ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
#     # 
#     # 
#     # 
#     # png(paste0(basepath, gene, ".png"), width = 300, height = 300)
#     # print(plot)
#     # dev.off()
#   }
# }
# 
# rsq_df_neurons <- rsq_df_neurons[order(-rsq_df_neurons$rsquare),]
# 
# topgenes_logodds_neurons <- c()
# for(i in 1 : 30) {
#   if (rsq_df_neurons$rsquare[i] >= 0.4) {
#     topgenes_logodds_neurons <- append(topgenes_logodds_neurons, rsq_df_neurons$gene[i])
#   }
# }
# 
# 
# training.data.topgenes.neurons <- training.data[, c(topgenes_logodds_neurons, "Main_cluster_name")]
# testing.data.topgenes.neurons <- testing.data[, c(topgenes_logodds_neurons, "Main_cluster_name")]
# 
# 
# training.data.topgenes.neurons$celltype <- 0
# training.data.topgenes.neurons$celltype[training.data.topgenes.neurons$Main_cluster_name == "Granule neurons"] <- 1
# training.data.topgenes.neurons$celltype[training.data.topgenes.neurons$Main_cluster_name == "Inhibitory interneurons"] <- 1
# 
# 
# testing.data.topgenes.neurons$celltype <- 0
# testing.data.topgenes.neurons$celltype[testing.data.topgenes.neurons$Main_cluster_name == "Granule neurons"] <- 1
# testing.data.topgenes.neurons$celltype[testing.data.topgenes.neurons$Main_cluster_name == "Inhibitory interneurons"] <- 1
# 
# 
# training.data.topgenes.neurons$Main_cluster_name <- NULL
# testing.data.topgenes.neurons$Main_cluster_name <- NULL
# 
# 
# 
# 
# #training glm with 0.75 GE data
# model_lm_neurons <- glm(celltype ~. , data = training.data.topgenes.neurons, family = "binomial")
# 
# #testing GLM with 0.25 GE data
# predictions.neurons <- model_lm_neurons %>% predict(testing.data.topgenes.neurons)
# 
# #determining the model performance
# performance.neurons <- data.frame(RMSE = RMSE(predictions.neurons, testing.data.topgenes.neurons$celltype),
#                                   R2 = R2(predictions.neurons, testing.data.topgenes.neurons$celltype))
# 
# 
# #view model summmary and performance
# summary(model_lm_neurons)
# performance.neurons
# 
# #vif to determine multicollinearity
# car::vif(model_lm_neurons)
# 
# #stepwise AIC 
# #stepwise regression --forward
# library(MASS)
# 
# 
# step.model.forward.neurons <- stepAIC(model_lm_neurons, direction = "forward", trace = 5)
# summary(step.model.forward.neurons)
# 
# step.model.backward.neurons <- stepAIC(model_lm_neurons, direction = "backward", trace = 5)
# summary(step.model.backward.neurons)
# 
# step.model.both.neurons <- stepAIC(model_lm_neurons, direction = "both", trace = 5)
# summary(step.model.both.neurons)


day110_gene_expression_data$celltype <- 0
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Granule neurons"] <- 1
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Inhibitory interneurons"] <- 1

#day 110 predictions neuron model
# Making predictions
probabilities.neurons <- model_lm_neurons %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes.neurons <- ifelse(probabilities > 0.6, 1, 0)


# determining accuracy
table(predicted.classes.neurons == day110_gene_expression_data$celltype)
confusionMatrix(table(predicted.classes.neurons, day110_gene_expression_data$celltype))














#testing with day110 astrocyte model


day110_gene_expression_data <- as.data.frame(seuobj110@assays$RNA$data)

# switching cols and rows by transposing
day110_gene_expression_data <- t(day110_gene_expression_data)


# make data remain as data.frame
day110_gene_expression_data <- data.frame(day110_gene_expression_data)


# round vals to 2 spots
day110_gene_expression_data <- (round(day110_gene_expression_data, 2))


# add cell id col to metadata of seuobj
seuobj110$cellid <- rownames(seuobj110@meta.data)


# making a metadata dataframe of seuobj with main cluster and cell id
metadata_df <- seuobj110@meta.data[, c("Main_cluster_name", "cellid")]


# adding cell ids to metadata
day110_gene_expression_data$cellid <- rownames(day110_gene_expression_data)


# combining data and metadata
day110_gene_expression_data <- cbind(day110_gene_expression_data, metadata_df)


# removing the cell id from the expression dataframe
day110_gene_expression_data$cellid <- NULL


day110_gene_expression_data$celltype <- 0
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Astrocytes"] <- 1

# Making predictions
probabilities <- model_lm %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes <- ifelse(probabilities > 0.6, 1, 0)


# determining accuracy
table(predicted.classes == day110_gene_expression_data$celltype)
confusionMatrix(table(predicted.classes, day110_gene_expression_data$celltype))