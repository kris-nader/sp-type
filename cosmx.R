library(Seurat)

nano.obj <- LoadNanostring(data.dir = "data/nsclc_cosmx/Lung5_Rep1", fov = "lung5.rep1")

#add in precomputed Azimuth annotations
azimuth.data <- readRDS("data/nsclc_cosmx/nanostring_data.Rds")
nano.obj <- AddMetaData(nano.obj, metadata = azimuth.data$annotations)
nano.obj[["proj.umap"]] <- azimuth.data$umap
Idents(nano.obj) <- nano.obj$predicted.annotation.l1

# set to avoid error exceeding max allowed size of globals
options(future.globals.maxSize = 8000 * 1024^2)
nano.obj <- SCTransform(nano.obj, assay = "Nanostring", clip.range = c(-10, 10), verbose = FALSE)

# text display of annotations and prediction scores
head(slot(object = nano.obj, name = "meta.data")[2:5])

#load cosmx data clustered wth sctype and SCTypeDB


DimPlot(nano.obj, raster = FALSE)
ImageDimPlot(nano.obj, fov = "lung5.rep1", axes = FALSE, cols = "glasbey", group.by ='customclassif')

merged <- merge(data, nano.obj)
VariableFeatures(merged) <- c(VariableFeatures(data), VariableFeatures(nano.obj))
merged <- RunPCA(merged)
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged)
merged <- RunUMAP(merged, dims = 1:30)

#DimPlot(merged, reduction = "umap", group.by = "technique", raster = FALSE)
DimPlot(merged, reduction = "umap", group.by = "customclassif", raster = FALSE)
