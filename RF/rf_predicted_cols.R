library(ROCR)
library(car)

seuobj110 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj110.RDS")
astrocyte_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/astrocyte_rf.RDS")
olig_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Oligodendrocyte_rf.RDS")
ii_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Inhibitory_interneuronrf.RDS")
microglia_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Microglia_rf.RDS")
gn_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Granule_neuronrf.RDS")
pn_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Purkinje_neuron_rf.RDS")
ubc_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Unipolar_brush_cell_rf.RDS")
vec_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/Vascular_endothelial_cell_rf.RDS")
slc_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/RF/SLC24A4_PEX5L_positive_cell_rf.RDS")

models <- list("Astrocytes" = astrocyte_model, 
               "Inhibitory interneurons" = ii_model, 
               "Unipolar brush cells" = ubc_model, 
               "Oligodendrocytes" = olig_model,
               "Granule neurons" = gn_model,
               "Vascular endothelial cells" = vec_model,
               "Microglia" = microglia_model,
               "Purkinje neurons" = pn_model,
               "SLC24A4_PEX5L positive cells" = slc_model)

#day 110 data frame

#testing with day110 astrocyte model

day110_gene_expression_data <- as.data.frame(seuobj110@assays$RNA$data)
day110_gene_expression_data <- mutate_all(day110_gene_expression_data, function(x) as.numeric(as.character(x)))
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





#adjusted probability calc 




# 
# probabilities <- model_lm %>% predict(day110_gene_expression_data, type = "response") 
# 
original_odds <- 72266/217612
# 
# undersample_odds <- sum(training.data$Main_cluster_name == "Astrocytes")/sum(training.data$Main_cluster_name != "Astrocytes")
# 
# scoring_odds <- probabilities/(1-probabilities)
# 
# adjusted_odds <- scoring_odds * (original_odds/undersample_odds)
# 
# adjusted_probability = 1/(1+(1/adjusted_odds))

seuobj_115_125 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/seuobj_115_125.RDS")

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








for(i in names(models)){
  
  
  
  probabilities <- models[[i]] %>% predict(day110_gene_expression_data, type = "response")
  undersample_odds <- sum(training.data$Main_cluster_name == i)/sum(training.data$Main_cluster_name != i)
  scoring_odds <- probabilities/(1-probabilities)
  adjusted_odds <- scoring_odds * (original_odds/undersample_odds)
  adjusted_probability = 1/(1+(1/adjusted_odds))
  colname <- paste(i, "probability")
  seuobj110@meta.data[, colname] <- 1 - adjusted_probability
  
  
  
}


df <- seuobj110@meta.data[,grep(" probability", colnames(seuobj110@meta.data))]
df <- colnames(df)[apply(df,1,which.max)]

seuobj110@meta.data$max_probability <- df 


df3 <- seuobj110@meta.data[,grep("Day89_Topic_", colnames(seuobj110@meta.data))]
df3 <- colnames(df3)[apply(df3,1,which.max)]

seuobj110@meta.data$max_topic_day89 <- df3


df4 <- seuobj110@meta.data[,grep("Day94_Topic_", colnames(seuobj110@meta.data))]
df4 <- colnames(df4)[apply(df4,1,which.max)]

seuobj110@meta.data$max_topic_day94 <- df4


df5 <- seuobj110@meta.data[,grep("Topic__", colnames(seuobj110@meta.data))]
df5 <- colnames(df5)[apply(df5,1,which.max)]

seuobj110@meta.data$max_topic_day110 <- df5


df6 <- seuobj110@meta.data[,grep("Day115_Topic_", colnames(seuobj110@meta.data))]
df6 <- colnames(df6)[apply(df6,1,which.max)]

seuobj110@meta.data$max_topic_day115 <- df6


df7 <- seuobj110@meta.data[,grep("Day125_Topic_", colnames(seuobj110@meta.data))]
df7 <- colnames(df7)[apply(df7,1,which.max)]

seuobj110@meta.data$max_topic_day125 <- df7


plotdf <- table(seuobj110$Main_cluster_name, seuobj110$max_probability)

plotdf <- as.data.frame(plotdf)

ggplot(plotdf, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs predicted cluster name in Day 110")



plotdf <- table(seuobj110$Main_cluster_name, seuobj110$max_topic_day89)

plotdf <- as.data.frame(plotdf)
plotdf$Var2 <- gsub("Day89_", "", plotdf$Var2)
plotdf$Var2 <- factor(plotdf$Var2, levels = 
                        mixedsort(unique(plotdf$Var2)))

ggplot(plotdf, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 89 topics in day 110 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))




plotdf2 <- table(seuobj110$Main_cluster_name, seuobj110$max_topic_day94)

plotdf2 <- as.data.frame(plotdf2)
plotdf2$Var2 <- gsub("Day94_", "", plotdf2$Var2)
plotdf2$Var2 <- factor(plotdf2$Var2, levels = 
                         mixedsort(unique(plotdf2$Var2)))

ggplot(plotdf2, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 94 topics in day 110 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))


plotdf3 <- table(seuobj110$Main_cluster_name, seuobj110$max_topic_day115)

plotdf3 <- as.data.frame(plotdf3)
plotdf3$Var2 <- gsub("Day115_", "", plotdf3$Var2)
plotdf3$Var2 <- factor(plotdf3$Var2, levels = 
                         mixedsort(unique(plotdf3$Var2)))

ggplot(plotdf3, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 115 topics in day 110 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))



plotdf4 <- table(seuobj110$Main_cluster_name, seuobj110$max_topic_day125)

plotdf4 <- as.data.frame(plotdf4)
plotdf4$Var2 <- gsub("Day125_", "", plotdf4$Var2)
plotdf4$Var2 <- factor(plotdf4$Var2, levels = 
                         mixedsort(unique(plotdf4$Var2)))

ggplot(plotdf4, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 125 topics in day 110 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))