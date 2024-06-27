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

cds <- cluster_cells(cds, resolution=1e-5)
plot_cells(cds)

# visualize partition
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
plot_cells(cds, color_cells_by="cao_cell_type")
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE)

# find marker genes expressed by each cluster
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

# plot expression and partitions
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)

# ------------------------------------------------------------------------------
# It's often informative to look at more than one marker, which you can do just by changing the first argument to top_n():
# top_specific_markers <- marker_test_res %>%
#                             filter(fraction_expressing >= 0.10) %>%
#                             group_by(cell_group) %>%
#                             top_n(3, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="cluster_row_col",
#                     max.size=3)

# annotate cells according to type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))

  # remap each cluster to a different cell type
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="Body wall muscle",
                                                 "2"="Germline",
                                                 "3"="Motor neurons",
                                                 "4"="Seam cells",
                                                 "5"="Sex myoblasts",
                                                 "6"="Socket cells",
                                                 "7"="Marginal_cell",
                                                 "8"="Coelomocyte",
                                                 "9"="Am/PH sheath cells",
                                                 "10"="Ciliated neurons",
                                                 "11"="Intestinal/rectal muscle",
                                                 "12"="Excretory gland",
                                                 "13"="Chemosensory neurons",
                                                 "14"="Interneurons",
                                                 "15"="Unclassified eurons",
                                                 "16"="Ciliated neurons",
                                                 "17"="Pharyngeal gland cells",
                                                 "18"="Unclassified neurons",
                                                 "19"="Chemosensory neurons",
                                                 "20"="Ciliated neurons",
                                                 "21"="Ciliated neurons",
                                                 "22"="Inner labial neuron",
                                                 "23"="Ciliated neurons",
                                                 "24"="Ciliated neurons",
                                                 "25"="Ciliated neurons",
                                                 "26"="Hypodermal cells",
                                                 "27"="Mesodermal cells",
                                                 "28"="Motor neurons",
                                                 "29"="Pharyngeal gland cells",
                                                 "30"="Ciliated neurons",
                                                 "31"="Excretory cells",
                                                 "32"="Amphid neuron",
                                                 "33"="Pharyngeal muscle")
plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type")

# look at one partition isolated
cds_subset <- choose_cells(cds)


