library(Seurat)

data <- Read10X("data/cancer/brca_xenvis/sc/ftp", gene.column = 2)
data <- CreateSeuratObject(data)
data <- SCTransform(data, ncells = 3000, conserve.memory = TRUE)
data <- RunPCA(data, assay = "SCT", verbose = TRUE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.8) #0.8 is the default resolution
data <- RunUMAP(data, reduction = "pca", dim = 1:30)

save(data, file = "brca_xenvis_sc.rds")

source("sptype_KMN.R")
load("rds/brca_xenvis_sc.rds")
load("rds/brca_ref_markers_minor_normalized.rds")
marker_per_cluster <- NULL
tissue <- "Breast"
retrieve_markers(data,
                 ref_markers,
                 tissue = tissue,
                 marker_per_cluster = marker_per_cluster)
data <- run_sctype(seurat_object = data,
           known_tissue_type = "Breast",
           custom_marker_file = "temp/sctypeDB_.xlsx",
           slot = "SCT")
DimPlot(data, group.by = "sctype_classification")
