library(Seurat)

load("brca.rds")
load("brca_ref.rds")

ref <- subset(ref, subtype == "HER2+")
ref <- SCTransform(ref, ncells = 3000)
ref <- RunPCA(ref)
ref <- RunUMAP(ref, dims = 1:30)

anchors <- FindTransferAnchors(reference = ref, query = data, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype_minor, prediction.assay = TRUE,
                                  weight.reduction = data[["pca"]], dims = 1:30)
data[["predictions"]] <- predictions.assay
DefaultAssay(data) <- "predictions"
data <- FindSpatiallyVariableFeatures(data, assay = "predictions", selection.method = "moransi",
                                        features = rownames(data), r.metric = 5, slot = "data")
top.clusters = rownames(dplyr::slice_min(data[["predictions"]]@meta.features,
                                         moransi.spatially.variable.rank, n = 15))
SpatialPlot(object = data, features = top.clusters[1:3])
