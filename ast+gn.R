seuobj110$predicted.classes_Astrocytes <- predicted.classes_astrocytes
seuobj110$predicted.classes_Granule_neurons <- predicted.classes_Granule_neurons

df <- seuobj110@meta.data[,grep("predicted.classes_1", colnames(seuobj110@meta.data))]
df <- colnames(df)[apply(df,1,which.max)]

seuobj110$max_probability <- df
plot_df <- table(seuobj110$Main_cluster_name, seuobj110$max_probability)

library(ggplot2)

# value:table
library(ROCR)
library(car)
seuobj110 <- readRDS("/Users/saanviiyer/PRINT OPTICS/STUDY/SCRI_r/fetal_cerebellar_scData/seuobj110.RDS")
astrocyte_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/astrocyte_model.RDS")
gn_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/0.4+/gn_model.RDS")

models <- list("Astrocytes" = astrocyte_model, 
               "Granule neurons" = gn_model)

models_str <- c("astrocyte", "GN")
cell_types <- c("Astrocytes", "Granule neurons")


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

original_odds <- 72266/217612

for(i in names(models)){
  
  probabilities <- models[[i]] %>% predict(day110_gene_expression_data, type = "response")
  undersample_odds <- sum(training.data$Main_cluster_name == i)/sum(training.data$Main_cluster_name != i)
  scoring_odds <- probabilities/(1-probabilities)
  adjusted_odds <- scoring_odds * (original_odds/undersample_odds)
  adjusted_probability = 1/(1+(1/adjusted_odds))
  colname <- paste(i, "probability")
  seuobj110@meta.data[, colname] <- adjusted_probability
  
}

df <- seuobj110@meta.data[,grep(" probability", colnames(seuobj110@meta.data))]
df <- colnames(df)[apply(df,1,which.max)]

seuobj110@meta.data$max_probability <- df 


plotdf <- table(seuobj110$Main_cluster_name, seuobj110$max_probability)

plotdf <- as.data.frame(plotdf)

ggplot(plotdf, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs predicted cluster name in Day 110")
