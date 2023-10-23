library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)

#read reference sc file
#ref from https://www.nature.com/articles/s41588-021-00911-1
ref <- Read10X("brca_ref", gene.column = 1)
#read metadata
metadata <- read.csv("brca_ref/metadata.csv")
cell_types <- metadata$celltype_minor
names(cell_types) <- metadata$X
cell_types <- as.factor(cell_types)
#create the Reference object
ref <- Reference(ref, cell_types)

#load spatial data
data_dir <- "data/brca"
filename <- "filtered_feature_bc_matrix.h5"
data <- Load10X_Spatial(data.dir = data_dir, filename = filename)
coords <- GetTissueCoordinates(data)
counts <- data@assays$Spatial@counts
#create SpatialRNA object
puck <- SpatialRNA(coords, counts)

rctd <- create.RCTD(puck, ref)
rctd <- run.RCTD(rctd, doublet_mode = 'full')

save(rctd, file = "brca_rctd_minor.rds")

load("brca_rctd_minor.rds")
results <- rctd@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- rctd@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- rctd@spatialRNA
resultsdir <- "figures/brca/rctd/minor"
# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

weight_matrix <- as.matrix(norm_weights)
#for each barcode get the cell type with highest weight
rctd_annots <- colnames(weight_matrix)[max.col(weight_matrix,ties.method="random")]
#load data and add annotation
data@meta.data$rctd_annots <- rctd_annots
#plot

SpatialDimPlot(data, group.by = 'rctd_annots',
                     pt.size.factor = 2.4,
                     alpha = c(1, 1), label = TRUE, label.size = 2) +
  guides(colour = guide_legend(override.aes = list(size=10)))
