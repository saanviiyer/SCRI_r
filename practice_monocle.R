# load libraries

library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)

# Load the data
expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))

# Make the CDS object | CDS = cell data set
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# pre-process the data | NORMALIZE DATA
cds <- preprocess_cds(cds, num_dim = 100, method = "PCA")

plot_pc_variance_explained(cds) # elbow plot

# reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds)

plot_cells(cds, color_cells_by="cao_cell_type")
# plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))

# check for and remove batch effects -> these are systematic differences in the transcriptome of cells measured in different experimental batches
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE) # plate annotation specifies which sci-RNA-seq plate each cell originated from

# remove batch effect through align_cds() function
cds <- align_cds(cds, num_dim = 100, alignment_group = "plate")
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

# visualize partition
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
plot_cells(cds, color_cells_by="cao_cell_type")

marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)



# plot expression and partitions
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)