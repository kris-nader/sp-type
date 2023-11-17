library(Seurat)

load("rds/brca_xenvis_visium.rds")

data <- SCTransform(data, assay="Spatial")
data <- RunPCA(data, assay = "SCT", verbose = TRUE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.8) #0.8 is the default resolution
data <- RunUMAP(data, reduction = "pca", dim = 1:30)

load("rds/brca_ref2.rds")

ref <- SCTransform(ref, ncells = 3000, conserve.memory = TRUE)
ref <- RunPCA(ref)
ref <- RunUMAP(ref, dims = 1:30)

save(ref, file = "rds/brca_ref2_normalized.rds")
load("rds/brca_ref2_normalized.rds")
seurat_time_3 <- system.time({
anchors <- FindTransferAnchors(reference = ref, query = data, normalization.method = "SCT", query.assay = "SCT", reference.assay = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE,
                                  weight.reduction = data[["pca"]], dims = 1:30)
})
data[["predictions"]] <- predictions.assay
DefaultAssay(data) <- "predictions"
data <- FindSpatiallyVariableFeatures(data, assay = "predictions", selection.method = "moransi",
                                        features = rownames(data), r.metric = 5, slot = "data")
top.clusters = rownames(dplyr::slice_min(data[["predictions"]]@meta.features,
                                         moransi.spatially.variable.rank, n = 15))
SpatialPlot(object = data, features = top.clusters[1:3])

weight_matrix <- data@assays$predictions@data %>% t()
weight_matrix <- weight_matrix %>% as.data.frame()
weight_matrix$max <- NULL
weight_matrix <- weight_matrix %>% as.matrix()
#for each barcode get the cell type with highest weight
data@meta.data$seurat_annots <- colnames(weight_matrix)[max.col(weight_matrix,ties.method="random")]

SpatialDimPlot(data, group.by = "seurat_annots", pt.size.factor = 2.4)
save(data, file = "rds/brca_xenvis_xenium.rds")
