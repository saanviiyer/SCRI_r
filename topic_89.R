library(dplyr)
seuobj89 <- readRDS("/Users/saanviiyer/Desktop/fetal_cerebellar_scData/seuobj89.RDS")
basepath <- "/Users/saanviiyer/Documents/GitHub/SCRI_r/topics"


df <- as.data.frame(seuobj89@meta.data)
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
  
  saveRDS(final_df, file = file.path(basepath, paste0(day, ".RDS")))
  
}

# -----------------------
# Add Day89 from Topic__ columns
# -----------------------
if (any(grepl("^Topic__", colnames(df)))) {
  day89_df <- df %>% select(starts_with("Topic__"))
  
  top_info <- t(apply(day89_df, 1, function(x) {
    ord <- order(x, decreasing = TRUE)[1:3]
    c(names(x)[ord], x[ord])
  }))
  
  top_df <- as.data.frame(top_info, stringsAsFactors = FALSE)
  colnames(top_df) <- c("Top1_name", "Top2_name", "Top3_name",
                        "Top1_value", "Top2_value", "Top3_value")
  
  day_dfs[["Day89"]] <- cbind(day89_df, top_df)
  
  saveRDS(final_df, file = file.path(basepath, "Day89.RDS"))
  
}

# Check what’s in the list
names(day_dfs)

# Load Day89 data with top topic columns
day89_topics <- day_dfs[["Day89"]]

# Match order of rows between Seurat metadata and Day89 dataframe
day89_topics <- day89_topics[match(rownames(seuobj89@meta.data), rownames(day89_topics)), ]

# Add top topic columns to Seurat metadata
seuobj89@meta.data$Top1_name  <- day89_topics$Top1_name
seuobj89@meta.data$Top2_name  <- day89_topics$Top2_name
seuobj89@meta.data$Top3_name  <- day89_topics$Top3_name
seuobj89@meta.data$Top1_value <- as.numeric(day89_topics$Top1_value)
seuobj89@meta.data$Top2_value <- as.numeric(day89_topics$Top2_value)
seuobj89@meta.data$Top3_value <- as.numeric(day89_topics$Top3_value)

# Subset to "Granule Neurons"
granule_obj <- subset(seuobj89, subset = Main_cluster_name == "Granule neurons")
# Create DimPlots without legends
p1 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top1_name") + NoLegend()
p2 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top2_name") + NoLegend()
p3 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top3_name") + NoLegend()

# Display plots
print(p1)
print(p2)
print(p3)

# Save plots
ggsave(file.path(basepath, "GranuleNeurons_Day89_Top1.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(basepath, "GranuleNeurons_Day89_Top2.png"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(basepath, "GranuleNeurons_Day89_Top3.png"), p3, width = 6, height = 5, dpi = 300)