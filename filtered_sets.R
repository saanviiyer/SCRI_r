seuobj110 <- subset(seuobj110, subset = Main_cluster_name == "Astrocytes" | Main_cluster_name == "Granule neurons")
astrocyte_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/astrocyte_model.RDS")
gn_model <- readRDS("/Users/saanviiyer/Documents/GitHub/SCRI_r/0.4+/gn_model.RDS")

probabilities <- astrocyte_model %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes_1ast <- ifelse(probabilities > 0.2, 1, 0)


probabilities <- gn_model %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes_1gn <- ifelse(probabilities > 0.2, 1, 0)


probabilities <- oligodendrocyte_model %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes_1olig <- ifelse(probabilities > 0.2, 1, 0)
