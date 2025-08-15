library(dplyr)
seuobj110 <- readRDS("/Users/saanviiyer/Desktop/fetal_cerebellar_scData/seuobj110.RDS")
basepath <- "/Users/saanviiyer/Documents/GitHub/SCRI_r/topics"


df <- as.data.frame(seuobj110@meta.data)
head(df)

# Extract the unique days from column names
days <- unique(sub("_.*", "", colnames(df)))

# Create a list to store each day's dataframe
day_dfs <- list()

# Extract the unique days from column names (only those starting with "Day")
days <- unique(sub("_.*", "", grep("^Day", colnames(df), value = TRUE)))

# Process normal DayXXX columns
for (day in days) {
  # Subset columns for this day
  day_df <- df %>% select(starts_with(day))
  
  # Get top 3 values and their names for each row
  top_info <- t(apply(day_df, 1, function(x) {
    ord <- order(x, decreasing = TRUE)[1:3]
    c(names(x)[ord], x[ord])
  }))
  
  # Convert to dataframe and add column names
  top_df <- as.data.frame(top_info, stringsAsFactors = FALSE)
  colnames(top_df) <- c("Top1_name", "Top2_name", "Top3_name",
                        "Top1_value", "Top2_value", "Top3_value")
  
  # Combine with original day dataframe
  final_df <- cbind(day_df, top_df)
  
  # Store in list
  day_dfs[[day]] <- final_df
  
  saveRDS(final_df, file = file.path(basepath, paste0(day, "_110.RDS")))
  
}

# -----------------------
# Add Day110 from Topic__ columns
# -----------------------
if (any(grepl("^Topic__", colnames(df)))) {
  day110_df <- df %>% select(starts_with("Topic__"))
  
  top_info <- t(apply(day110_df, 1, function(x) {
    ord <- order(x, decreasing = TRUE)[1:3]
    c(names(x)[ord], x[ord])
  }))
  
  top_df <- as.data.frame(top_info, stringsAsFactors = FALSE)
  colnames(top_df) <- c("Top1_name", "Top2_name", "Top3_name",
                        "Top1_value", "Top2_value", "Top3_value")
  
  day_dfs[["Day110"]] <- cbind(day110_df, top_df)
  
  saveRDS(final_df, file = file.path(basepath, "110_Day110.RDS"))
  
}

# Check what’s in the list
names(day_dfs)

# Load Day110 data with top topic columns
day110_topics <- day_dfs[["Day110"]]

# Match order of rows between Seurat metadata and Day110 dataframe
day110_topics <- day110_topics[match(rownames(seuobj110@meta.data), rownames(day110_topics)), ]

# Add top topic columns to Seurat metadata
seuobj110@meta.data$Top1_name  <- day110_topics$Top1_name
seuobj110@meta.data$Top2_name  <- day110_topics$Top2_name
seuobj110@meta.data$Top3_name  <- day110_topics$Top3_name
seuobj110@meta.data$Top1_value <- as.numeric(day110_topics$Top1_value)
seuobj110@meta.data$Top2_value <- as.numeric(day110_topics$Top2_value)
seuobj110@meta.data$Top3_value <- as.numeric(day110_topics$Top3_value)

# Subset to "Granule Neurons"
granule_obj <- subset(seuobj110, subset = Main_cluster_name == "Granule neurons")
# Create DimPlots without legends
p1 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top1_name") + NoLegend()
p2 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top2_name") + NoLegend()
p3 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top3_name") + NoLegend()

# Display plots
print(p1)
print(p2)
print(p3)

# Save plots
ggsave(file.path(basepath, "GranuleNeurons_Day110_Top1.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(basepath, "GranuleNeurons_Day110_Top2.png"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(basepath, "GranuleNeurons_Day110_Top3.png"), p3, width = 6, height = 5, dpi = 300)