library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(dplyr)
#read reference sc file
#ref from https://www.nature.com/articles/s41588-021-00911-1
ref <- Read10X("brca_ref", gene.column = 1)
#read metadata
metadata <- read.csv("brca_ref/metadata.csv")
cell_types <- metadata$celltype_minor
names(cell_types) <- metadata$X

#for brca xenvis
load("rds/brca_xenvis_visium.rds")
ref <- Read10X("data/cancer/brca_xenvis/sc/ftp", gene.column = 2)
metadata <- read.xlsx("data/cancer/brca_xenvis/Cell_Barcode_Type_Matrices.xlsx")
cell_types <- metadata$Annotation
names(cell_types) <- metadata$Barcode
cell_types <- as.factor(cell_types)
rctd_time_3 <- system.time({
#create the Reference object
ref <- Reference(ref, cell_types)

#load spatial data (raw)
#data_dir <- "data/cancer/brca_xenvis/visium"
#filename <- "filtered_feature_bc_matrix.h5"
#data <- Load10X_Spatial(data.dir = data_dir, filename = filename)
#load processed data

coords <- GetTissueCoordinates(data)
coords$cell <- NULL
counts <- data@assays$Spatial@counts #xenium or spatial

#create SpatialRNA object
puck <- SpatialRNA(coords, counts)

rctd <- create.RCTD(puck, ref)
rctd <- run.RCTD(rctd, doublet_mode = 'full')
})
save(rctd, file = "rds/brca_xenvis_xenium_rctd.rds")

load("rds/brca_xenvis_visium_rctd.rds")
results <- rctd@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- rctd@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- rctd@spatialRNA
resultsdir <- "figures/brca_xenvis/xenium/rctd"
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
rctd_annots <- data.frame(X = norm_weights@Dimnames[[1]],
                          rctd_annots = colnames(weight_matrix)[max.col(weight_matrix,ties.method="random")])
barcodes <- data.frame(X = data@assays$Xenium@counts@Dimnames[[2]]) #Xenium or Spatial
annots <- merge(rctd_annots, barcodes, all.y = TRUE, by = "X")
data@meta.data$rctd_annots <- annots$rctd_annots

save(data, file = "rds/brca_xenvis_xenium.rds")
#plot

SpatialDimPlot(data, group.by = 'rctd_annots',
                     pt.size.factor = 2.4,
                     alpha = c(1, 1))
