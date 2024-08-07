# load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(ggplot2)

# install.packages("ggpmisc")
# install.packages("broom")
library(broom)
library(ggpmisc)
library(nlme)

basepath <- "/Users/saanviiyer/Documents/GitHub/SCRI_r/0.4+/"

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
rsq_df

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
  if(nrow(count_df) >= 10){
    
    
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
list_genes <- c()

# for(gene in rsq_df$gene[1:30]){
for(i in 1 : 2000) {

  
  if (rsq_df$rsquare[i] >= 0.4) {
    print(paste0("COUNT",i))
    gene <- rsq_df$gene[i]
    message(gene)
    
    count <- unique(day115_125_gene_expression_data[,gene])
    
    count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
    
    message('Count df created')
    
    for(n in unique(count_df$count)){
      print(n)
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
      list_genes = append(gene, list_genes)
      # plot <- ggplot(count_df, aes(x = count, y = logodds)) +
      #   geom_point(stat = "identity") + geom_smooth(method = "lm") +
      #   ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
      # 
      # 
      # 
      # png(paste0(basepath, gene, "TOP30",".png"), width = 300, height = 300)
      # print(plot)
      # dev.off()
    }
    message("loop completed")
  }
}

list_genes
# 
# # creating top_genes_logodds
top_genes_logodds <- c()

# populating top_genes_logodds with top genes
for(i in 1 : 30) {

  if (rsq_df$rsquare[i] >= 0.2) {
    top_genes_logodds <- append(top_genes_logodds, rsq_df$gene[i])
  }
}

day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)

# split 75% training data
training.data <- day115_125_gene_expression_data %>%
  group_by(Main_cluster_name) %>%
  filter(row_number() <= 0.75 * n())

# all data not in training wil lbe used for testing
testing.data <- day115_125_gene_expression_data[!day115_125_gene_expression_data$cellid %in% training.data$cellid,]

rownames(training.data) <- training.data$cellid
rownames(testing.data) <- testing.data$cellid

training.data$cellid <- NULL
testing.data$cellid <- NULL

# # creating top genes from each
training.data.topgenes <- training.data[,c(top_genes_logodds, "Main_cluster_name")]
testing.data.topgenes <- testing.data[,c(top_genes_logodds, "Main_cluster_name")]

#getting the top genes from rsq_)df
topgenes_logodds <- c()
for(i in 1 : 30) {
  if (rsq_df$rsquare[i] >= 0.4) {
    topgenes_logodds <- append(topgenes_logodds, rsq_df$gene[i])
  }
}

newgenes <- c("S100B", "GFAP", "APOE", "AQP4")

for (gene in newgenes) {
  plot <- ggplot(count_df, aes(x = count, y = logodds)) +
    geom_point(stat = "identity") + geom_smooth(method = "lm") +
    ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
  
  
  
  png(paste0(basepath, gene, ".png"), width = 300, height = 300)
  print(plot)
  dev.off()
  
  topgenes_logodds <- append(topgenes_logodds, newgenes)
}

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

# NEURONAL MODEL


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
    
    plot <- ggplot(count_df, aes(x = count, y = logodds)) +
      geom_point(stat = "identity") + geom_smooth(method = "lm") +
      ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))



    png(paste0(basepath, gene, ".png"), width = 300, height = 300)
    print(plot)
    dev.off()
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

#getting the top genes from rsq_)df
topgenes_logodds <- c()
for(i in 1 : 2000) {
  if (rsq_df$rsquare[i] >= 0.4) {
    topgenes_logodds <- append(topgenes_logodds, rsq_df$gene[i])
  }
}

newgenes <- c("S100B", "GFAP", "APOE", "AQP4")

topgenes_logodds <- append(topgenes_logodds, newgenes[!newgenes %in% topgenes_logodds])

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
training.data.topgenes$celltype <- 1
training.data.topgenes$celltype[training.data.topgenes$Main_cluster_name == "Astrocytes"] <- 0

testing.data.topgenes$celltype <- 1
testing.data.topgenes$celltype[testing.data.topgenes$Main_cluster_name == "Astrocytes"] <- 0

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

# switching positives for model - predict astrocyte
day110_gene_expression_data$celltype <- 1
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Astrocytes"] <- 0

# Making predictions with step.model.both
probabilities <- step.model.both %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes <- ifelse(probabilities > 0.2, 1, 0)


# determining accuracy
table(predicted.classes == day110_gene_expression_data$celltype)
confusionMatrix(table(predicted.classes, day110_gene_expression_data$celltype))

#------------------------
# original odds = number of astrocytes / number of not-astrocytes in training data
orig_odds = (72266/217612)
# undersample odds <- number of astrocytes / number of not-astrocytes in undersampled data
undersampled_odds <- (sum(training.data$Main_cluster_name == "Astrocytes")/(sum(training.data$Main_cluster_name != "Astrocytes")))
# scoring odds = predicted probabilities / (1-predicted probabilities)
scoring_odds <- (probabilities) / (1-probabilities)
# adjusted odds = scoring odds * original odds / undersampled odds
adj_odds = scoring_odds * orig_odds / undersampled_odds
# adjusted probability = 1 / (1 + (1 / adjusted odds))
adj_prob = 1 / (1 + (1 / adj_odds))
# adjust our predictions with adjusted probabilities
predicted.classes <- ifelse(adj_prob > 0.2, 1, 0) # adj prob to change
conf_matrix <- confusionMatrix(table(predicted.classes, day110_gene_expression_data$celltype))

## calculating AUC / ROC
library(ROCR)
roc_testing_prediction <- prediction(adj_prob, day110_gene_expression_data$celltype)

# analyzing performance with true positive rate and false positive rate
roc_testing_performance <- performance(roc_testing_prediction, "tpr", "fpr")

# evaluating area under curve
auc_testing_performance <- performance(roc_testing_prediction, measure="auc")

saveRDS(step.model.both, file=paste0(basepath, "astrocyte_model.RDS"))

auc_testing_performance <- auc_testing_performance@y.values

# metrics for further analysis
tuning_metrics_df <- data.frame(cutoffs=NA, f_measure=NA, precision=NA, recall=NA)

TuningMetrics.fn <- function(beta, adjusted_probability, seuobj110){
  for(cutoff in seq(0.2, 0.9, 0.05)){ # switch 0.1 to 0.01
    predicted.classes <- ifelse(adjusted_probability > cutoff, 1, 0 )
    conf_matrix <- confusionMatrix(table(predicted.classes, seuobj110[["celltype"]]))
    
    TP <- conf_matrix$table[4]
    TN <- conf_matrix$table[1]
    FP <- conf_matrix$table[3]
    FN <- conf_matrix$table[2]
    
    precision <- TP/(TP + FP)
    recall <- TP/ (TP + FN)
    
    f_measure <- (1+beta**2) * ( (precision * recall) / ((beta**2 * precision) + recall) )
    
    
    tuning_metrics_df <- rbind(tuning_metrics_df, c(cutoff, f_measure, precision, recall))
    return(tuning_metrics_df)
  }
}

tuning_metrics_df_b01 <- TuningMetrics.fn(beta=0.1, adj_prob, day110_gene_expression_data)
df <- tuning_metrics_df_b01 %>% gather(metric, value, f_measure:recall)
ggplot(df, aes(x=cutoffs, y=value, color=metric)) + geom_point()

tuning_metrics_df_b01 <- TuningMetrics.fn(beta=0.1, adj_prob, day110_gene_expression_data)
tuning_metrics_df_b05 <- TuningMetrics.fn(beta=0.5, adj_prob, day110_gene_expression_data)
tuning_metrics_df_b1 <- TuningMetrics.fn(beta=1, adj_prob, day110_gene_expression_data)
tuning_metrics_df_b1.5 <- TuningMetrics.fn(beta=1.5, adj_prob, day110_gene_expression_data)
tuning_metrics_df_b2 <- TuningMetrics.fn(beta=2, adj_prob, day110_gene_expression_data)

seuobj110@meta.data$predicted_astrocytes_temp <- predicted.classes

# switching from wide to long dataframe
df <- tuning_metrics_df %>% gather(metric, value, f_measure:recall)
tuning_metrics_df_b05 <- tuning_metrics_df %>% gather(metric, value, f_measure:recall)
tuning_metrics_df_b1 <- tuning_metrics_df %>% gather(metric, value, f_measure:recall)
tuning_metrics_df_b2 <- tuning_metrics_df %>% gather(metric, value, f_measure:recall)

ggplot(tuning_metrics_df_b05, aes(x=cutoffs, y=value, color=metric)) + geom_point()
ggplot(tuning_metrics_df_b1, aes(x=cutoffs, y=value, color=metric)) + geom_point()
ggplot(tuning_metrics_df_b2, aes(x=cutoffs, y=value, color=metric)) + geom_point()

predicted.classes <- (adj_prob > 0.2)
seuobj110@meta.data$predicted_astrocytes_temp <- predicted.classes
DimPlot(seuobj110, reduction="umap", group.by = "predicted_astrocytes_temp", label = TRUE)

seuobj110@meta.data$predicted_astrocytes_probabilities <- predicted.classes
predicted.classes_astrocytes <- predicted.classes