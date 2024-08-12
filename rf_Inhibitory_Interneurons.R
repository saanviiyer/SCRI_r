library(ggplot2)
library(ggpmisc)
library(nlme)
library(broom)
library(dplyr)
library(car)
library(Seurat)
library(randomForest)
library(datasets)

seuobj110 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj110.RDS")


seuobj_115_125 <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/day115_125.RDS")

#creating data frame of gene expression data
day115_125_gene_expression_data <- as.data.frame(seuobj_115_125@assays$integrated$data)

# switching cols and rows by transposing
day115_125_gene_expression_data <- t(day115_125_gene_expression_data)


# make data remain as data.frame
day115_125_gene_expression_data <- data.frame(day115_125_gene_expression_data)


# round vals to 2 spots
day115_125_gene_expression_data <- (round(day115_125_gene_expression_data, 2))


# add cell id col to metadata of seuobj
seuobj_115_125$cellid <- rownames(seuobj_115_125@meta.data)


# making a metadata dataframe of seuobj with main cluster and cell id
metadata_df <- seuobj_115_125@meta.data[, c("Main_cluster_name", "cellid")]


# adding cell ids to metadata
day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)


# combining data and metadata
day115_125_gene_expression_data <- cbind(day115_125_gene_expression_data, metadata_df)


# removing the cell id from the expression dataframe
day115_125_gene_expression_data$cellid <- NULL


#create a data frame to store r^2 values for each gene
rsq_df <- data.frame(gene = colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"], rsquare = rep(NA, 2001))

rsq_df[]




#for loop to populate rsq_df
for(gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"]){
  #tells name of current gene
  message(gene)
  
  #filling count with gene expression values for all cells for this gene
  count <- unique(day115_125_gene_expression_data[,gene])
  
  #creating a dataframe with percent and logodds columns corresponding to probability that cell is an Inhibitory_interneuron depending on counts for that gene
  count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
  
  #displays that coutn df was created
  message('Count df created')
  
  #for loop iterating through each unique count value
  for(n in unique(count_df$count)){
    
    # count for all Inhibitory_interneurons
    n_Inhibitory_interneurons <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Inhibitory interneurons",gene]==n)
    
    # probability of a cell being an Inhibitory_interneuron
    # number of Inhibitory_interneurons w unique count value for spp1 / all cells
    percent <- n_Inhibitory_interneurons/ sum(day115_125_gene_expression_data[,gene]==n)
    
    # log odds
    logodds <- log(percent/(1-percent))
    
    # adding all variables to count_df
    count_df[count_df$count == n, "percent"] <- percent
    count_df[count_df$count == n, "logodds"] <- logodds
  }
  
  
  #removing all -Inf logodds
  count_df <- count_df[count_df$logodds != -Inf, ]
  #removing all Inf logodds
  count_df <- count_df[count_df$logodds != Inf, ]
  
  #only executing if there are more than 20 unique counts for the gene
  if(nrow(count_df) >= 10){
    
    #creating a linear model to model relation between count and logodds
    model <- lm(logodds ~ count, data = count_df)
    #getting summary of the model
    model_summary <- glance(model)
    #extracting r squared from the summary
    r_squared <- model_summary$r.squared
    
    #adding r62 value from the model into rsq_df
    rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared
  }
}




#ordering rsq_df by descending r^2 
rsq_df <- rsq_df[order(-rsq_df$rsquare),]




#Inhibitory_interneuron model
#getting the top genes from rsq_)df
topgenes_logodds <- c()
for(i in 1 : length(rsq_df)) {
  if (rsq_df$rsquare[i] >= 0.4) {
    topgenes_logodds <- append(topgenes_logodds, rsq_df$gene[i])
  }else{
    break
  }
}

Inhibitory_interneuron.markers <- c("S100B","GFAP", "APOE", "AQP4")
Inhibitory_interneuron.markers <- Inhibitory_interneuron.markers[!Inhibitory_interneuron.markers %in% topgenes_logodds]
topgenes_logodds <- append(topgenes_logodds, Inhibitory_interneuron.markers)



#load dplyr
library(dplyr)

#making the cell id columns the same as the rownames
day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)

#making the training dat groping by main cluste name and filtering out
training.data <- day115_125_gene_expression_data %>%
  group_by(Main_cluster_name) %>%
  filter(row_number() <= 0.75 * n())

#making testing data from everything that isnt training
testing.data <- day115_125_gene_expression_data[!day115_125_gene_expression_data$cellid %in% training.data$cellid, ]


rownames(training.data) <- training.data$cellid
rownames(testing.data) <- testing.data$cellid


training.data$cellid <- NULL
testing.data$cellid <- NULL


#creating top genes from each
training.data.topgenes <- training.data[, c(topgenes_logodds, "Main_cluster_name")]
testing.data.topgenes <- testing.data[, c(topgenes_logodds, "Main_cluster_name")]

#creating binary predcictor column (Inhibitory_interneuron = 1, other = 0)
training.data.topgenes$celltype <- 1
training.data.topgenes$celltype[training.data.topgenes$Main_cluster_name == "Inhibitory interneurons"] <- 0

testing.data.topgenes$celltype <- 1
testing.data.topgenes$celltype[testing.data.topgenes$Main_cluster_name == "Inhibitory interneurons"] <- 0




training.data.topgenes$Main_cluster_name <- NULL
testing.data.topgenes$Main_cluster_name <- NULL


#training glm with 0.75 GE data
model_lm <- randomForest(celltype ~. , data = training.data.topgenes, family = "binomial")

#testing GLM with 0.25 GE data
predictions <- model_lm %>% predict(testing.data.topgenes)

#determining the model performance
performance <- data.frame(RMSE = RMSE(predictions, testing.data.topgenes$celltype),
                          R2 = R2(predictions, testing.data.topgenes$celltype))


#view model summmary and performance
summary(model_lm)
performance




#stepwise AIC 
#stepwise regression --forward
library(MASS)


# step.model.forward <- stepAIC(model_lm, direction = "forward", trace = 5)
# summary(step.model.forward)
# 
# step.model.backward <- stepAIC(model_lm, direction = "backward", trace = 5)
# summary(step.model.backward)
# 
# step.model.both <- stepAIC(model_lm, direction = "both", trace = 5)
# summary(step.model.both)













#testing with day110 Inhibitory_interneuron model
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


day110_gene_expression_data$celltype <- 1
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Inhibitory interneurons"] <- 0

# Making predictions
probabilities <- model_lm %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes <- ifelse(probabilities > 0.6, 1, 0)


# determining accuracy
table(predicted.classes == day110_gene_expression_data$celltype)
confusionMatrix(table(predicted.classes, day110_gene_expression_data$celltype))


#original odds = number of Inhibitory_interneurons / number of not Inhibitory_interneurons in training data
#undersample odds <- number of Inhibitory_interneurons / number of not Inhibitory_interneurons in undersampled data
#scoring odds <- predicted probabilities / (1-predicted probabilties)
#adjusted odds = scoring odds * original odds/undersampled odds
#adjusted probability = 1 / (1+(1/adjusted odds))

original_odds <- 72266/217612

undersample_odds <- sum(training.data$Main_cluster_name == "Inhibitory interneurons")/sum(training.data$Main_cluster_name != "Inhibitory interneurons")

scoring_odds <- probabilities/(1-probabilities)

adjusted_odds <- scoring_odds * (original_odds/undersample_odds)

adjusted_probability = 1/(1+(1/adjusted_odds))







#adjust our predictions with adjusted probabilites
predicted.classes <- ifelse(adjusted_probability > 0.6, 1, 0) #adjusted prob to change
conf_matrix <- confusionMatrix(table(predicted.classes, day110_gene_expression_data$celltype))

#calculating AUC/ROC
roc_testing_prediction <- prediction(adjusted_probability, day110_gene_expression_data$celltype)
roc_testing_performance <- performance(roc_testing_prediction, "tpr", "fpr")
auc_testing_performance <- performance(roc_testing_prediction, measure = "auc")
auc_testing_performance <- auc_testing_performance@y.values



tuning_metrics_df <- data.frame(cutoffs = NA, f_measure = NA, precision = NA, recall = NA)

#beta value for calculating f measure
beta <- 0.5

#iterating through different cutoff values for the predicted classes and calculating metrics for each value
for(cutoff in seq(0.2, 1, 0.1)) { # switch 0.1 0.01
  #make the predicted classes
  predicted.classes <- ifelse(adjusted_probability > cutoff, 1, 0)
  #make a confusion matrix with the predicted classes and the actual day110 cell types
  conf_matrix <- confusionMatrix(table(predicted.classes, day110_gene_expression_data$celltype))
  
  #assigning values from confusion matrix
  TP <- conf_matrix$table[4]
  TN <- conf_matrix$table[1]
  FP <- conf_matrix$table[3]
  FN <- conf_matrix$table[2]
  
  
  #calculating precision
  precision <- TP/(TP + FP)
  #calculating recall
  recall <- TP/(TP + FN)
  #calculating f measure
  f_measure <- ((1+beta**2) * (precision * recall)) / ((beta**2 * precision) + recall)
  
  #adding calculated values for current cutoff to tuning metrics df
  tuning_metrics_df <- rbind(tuning_metrics_df, c(cutoff, f_measure, precision, recall))
  
}

#save the nde
saveRDS(model_lm, file = "/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Inhibitory_interneuronrf.RDS")

#setting tuning metrics to the one for the regular model
tuning_metrics_df_model_lm <- tuning_metrics_df

#switching from wide to long dataframe
df <- tuning_metrics_df_model_lm %>% gather(metric, value, f_measure:recall)

#make ggplot of the data
ggplot(df, aes(x = cutoffs, y = value, color = metric)) + geom_point()

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

tuning_metrics_df_b01 <- TuningMetrics.fn(beta=0.1, adjusted_probability, day110_gene_expression_data)
df <- tuning_metrics_df_b01 %>% gather(metric, value, f_measure:recall)
ggplot(df, aes(x=cutoffs, y=value, color=metric)) + geom_point()

tuning_metrics_df_b01 <- TuningMetrics.fn(beta=0.1, adjusted_probability, day110_gene_expression_data)
tuning_metrics_df_b05 <- TuningMetrics.fn(beta=0.5, adjusted_probability, day110_gene_expression_data)
tuning_metrics_df_b1 <- TuningMetrics.fn(beta=1, adjusted_probability, day110_gene_expression_data)
tuning_metrics_df_b1.5 <- TuningMetrics.fn(beta=1.5, adjusted_probability, day110_gene_expression_data)
tuning_metrics_df_b2 <- TuningMetrics.fn(beta=2, adjusted_probability, day110_gene_expression_data)

df_b05 <- tuning_metrics_df_b05 %>% gather(metric, value, f_measure:recall)
df_b1 <- tuning_metrics_df_b1 %>% gather(metric, value, f_measure:recall)
df_b2 <- tuning_metrics_df_b2 %>% gather(metric, value, f_measure:recall)

ggplot(df_b05, aes(x = cutoffs, y = value, color = metric)) + geom_point()
ggplot(df_b1, aes(x = cutoffs, y = value, color = metric)) + geom_point()
ggplot(df_b2, aes(x = cutoffs, y = value, color = metric)) + geom_point()




predicted.classes <- (adjusted_probability > 0.9)



# 
# for(gene in Inhibitory_interneuron.markers){
#   plot <- ggplot(count_df, aes(x = count, y = logodds)) +
#     geom_point(stat = "identity") + geom_smooth(method = "lm") +
#     ggtitle(paste0("Log odds for ", gene, "\n(Inhibitory_interneurons)")) + stat_poly_eq(use_label(c("eq", "R2")), formula = y~x)
#   
#   
#   
#   png(paste0(basepath, gene, ".png"), width = 300, height = 300)
#   print(plot)
#   dev.off()
# }


seuobj110$probabilities <- probabilities

seuobj110$predicted_Inhibitory_interneurons_temp <- predicted.classes

seuobj110$probabilities_neg <- 1-probabilities
DimPlot(seuobj110, group.by = "predicted_Inhibitory_interneurons_temp")
FeaturePlot(seuobj110, features = "probabilities_neg", min.cutoff = "q1")