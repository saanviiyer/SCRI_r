#load libraries
library(monocle3)
library(ggplot2)
library(dplyr)

# Load the data
expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))

#CDS object
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

#pre-processing
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

#dimensionality reduction pca
cds <- reduce_dimension(cds, preprocess_method = 'PCA')
#plot cells
plot_cells(cds)
#plot manually annotated cell types
plot_cells(cds, color_cells_by="cao_cell_type")
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))
#dimensionality reduction tSNE
cds <- reduce_dimension(cds, reduction_method="tSNE")
#tSNE plot
plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type")
#color by plate (batch effects)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

#removing batch effects, aligns cells from different batches to allow for better comparison
cds <- align_cds(cds, num_dim = 100, alignment_group = "plate")
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

#clustering cells
cds <- cluster_cells(cds, resolution=1e-5)
plot_cells(cds)

#plot cells by partition
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

#add cell type annotations to partition plot
plot_cells(cds, color_cells_by="cao_cell_type")
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE)

#find marker genes expressed by each cluster
#tells us what identifies a cluster to be attributed to a given cell type
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)