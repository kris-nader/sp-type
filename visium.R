library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(HDF5Array)

#load data
data_dir <- "data/healthy/liver_crost_VISDP000066/350frompaper"
filename <- "filtered_feature_bc_matrix.h5"
data <- Load10X_Spatial(data.dir = data_dir, filename = filename)

#load data with no h5 file (brca wu dataset)
data_dir <- "data/human_intestine/A1"
mat_dir <- paste0(data_dir, "/raw_feature_bc_matrix")
im_dir <- paste0(data_dir, "/spatial")
mat <- Read10X(mat_dir, gene.column = 2)
data <- CreateSeuratObject(mat, assay = "Spatial")
im <- Read10X_Image(im_dir)
im <- im[Cells(x = data)]
DefaultAssay(object = im) <- 'Spatial'
data[['Slice1']] <- im
metadata <- read.csv(paste0(data_dir, "/metadata.csv"))
data@meta.data$subtype <- metadata$subtype
data@meta.data$patho_annot <- metadata$Classification
data@meta.data$patientid <- metadata$patientid

#quality control
#nFeature is the number of genes detected in each cell
#nCount is the total number of molecules detected within a cell
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
plot1 <- VlnPlot(data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                 pt.size = 0.1)
plot2 <- SpatialFeaturePlot(data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                            pt.size.factor = 2.4, alpha = c(0.3, 1))
#save quality control plots
output_folder <- "figures/human_intestine"
dir.create(output_folder, recursive = TRUE)
output_name <- "violin"
png(file = paste0("./", output_folder, "/", output_name, ".png"),
    width=1000, height=550)
print(plot1)
dev.off()
output_name <- "quality_control"
png(file = paste0("./", output_folder, "/", output_name, ".png"),
    width=1000, height=550)
print(plot2)
dev.off()
plot2
plot1
#filter
data <- data[, data$nFeature_Spatial > 300 & data$percent.mt < 25]

#normalization
data <- SCTransform(data, assay="Spatial")
#without using SCT:
#data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000) %>% ScaleData() %>% FindVariableFeatures()

#dim reduction and clustering
data <- RunPCA(data, assay = "SCT", verbose = TRUE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.8) #0.8 is the default resolution
data <- RunUMAP(data, reduction = "pca", dim = 1:30)

#save data as rds
save(data, file = "lymph_node_10x.rds")

#sptype
source("sptype_KMN.R")
data <- run_plot_sctype(data = data,
                        #db_ = "temp/sctypeDB_20.xlsx",
                        tissue = "Intestine",
                        output_folder = "figures/human_intestine",
                        output_name = "sptype",
                        saveRDS = FALSE,
                        pt.size.factor = 2)
save(data, file = "human_intestine.rds")
#plot clusters separately
SpatialDimPlot(data,
               cells.highlight = CellsByIdentities(object = data,
                                                   idents = c(1)),
               facet.highlight = TRUE, ncol = 6,
               pt.size.factor = 2.4)
#interactive linked plot 
LinkedDimPlot(data)

#markers for ident.1th cluster (if ident.2 is null, compares ident.1 to all other clusters)
markers <- FindMarkers(data, ident.1 = 1, ident.2 = NULL)
#plot top markers
SpatialFeaturePlot(object = data, features = rownames(markers)[1:3],
                   pt.size.factor = 2.4, alpha = c(0.1, 1), ncol = 3)

all_markers <- FindAllMarkers(data)
save(all_markers, file = "brca_markers.RData")

#get the marker with least p-val of each cluster
top_markers <- all_markers %>% arrange(p_val) %>% group_by(cluster) %>% slice(1)
row.names(top_markers) <- top_markers$gene
SpatialFeaturePlot(object = data, features = rownames(top_markers),
                   pt.size.factor = 2.4, alpha = c(0.1, 1), ncol = 5)

#get image in Seurat data (if image not available separately)
im <- GetImage(data, mode = "raster")
plot(im)