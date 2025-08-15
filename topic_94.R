library(dplyr)
seuobj94 <- readRDS("/Users/saanviiyer/Desktop/fetal_cerebellar_scData/seuobj94.RDS")
basepath <- "/Users/saanviiyer/Documents/GitHub/SCRI_r/topics"


df <- as.data.frame(seuobj94@meta.data)
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
  
  saveRDS(final_df, file = file.path(basepath, paste0(day, "_94.RDS")))
  
}

# -----------------------
# Add Day94 from Topic__ columns
# -----------------------
if (any(grepl("^Topic__", colnames(df)))) {
  day94_df <- df %>% select(starts_with("Topic__"))
  
  top_info <- t(apply(day94_df, 1, function(x) {
    ord <- order(x, decreasing = TRUE)[1:3]
    c(names(x)[ord], x[ord])
  }))
  
  top_df <- as.data.frame(top_info, stringsAsFactors = FALSE)
  colnames(top_df) <- c("Top1_name", "Top2_name", "Top3_name",
                        "Top1_value", "Top2_value", "Top3_value")
  
  day_dfs[["Day94"]] <- cbind(day94_df, top_df)
  
  saveRDS(final_df, file = file.path(basepath, "94_Day94.RDS"))
  
}

# Check what’s in the list
names(day_dfs)

# Load Day94 data with top topic columns
day94_topics <- day_dfs[["Day94"]]

# Match order of rows between Seurat metadata and Day94 dataframe
day94_topics <- day94_topics[match(rownames(seuobj94@meta.data), rownames(day94_topics)), ]

# Add top topic columns to Seurat metadata
seuobj94@meta.data$Top1_name  <- day94_topics$Top1_name
seuobj94@meta.data$Top2_name  <- day94_topics$Top2_name
seuobj94@meta.data$Top3_name  <- day94_topics$Top3_name
seuobj94@meta.data$Top1_value <- as.numeric(day94_topics$Top1_value)
seuobj94@meta.data$Top2_value <- as.numeric(day94_topics$Top2_value)
seuobj94@meta.data$Top3_value <- as.numeric(day94_topics$Top3_value)

# Subset to "Granule Neurons"
granule_obj <- subset(seuobj94, subset = Main_cluster_name == "Granule neurons")
# Create DimPlots without legends
p1 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top1_name") + NoLegend()
p2 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top2_name") + NoLegend()
p3 <- DimPlot(granule_obj, reduction = "umap", group.by = "Top3_name") + NoLegend()

# Display plots
print(p1)
print(p2)
print(p3)

# Save plots
ggsave(file.path(basepath, "GranuleNeurons_Day94_Top1.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(basepath, "GranuleNeurons_Day94_Top2.png"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(basepath, "GranuleNeurons_Day94_Top3.png"), p3, width = 6, height = 5, dpi = 300)

