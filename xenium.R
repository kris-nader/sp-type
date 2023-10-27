library(Seurat)

path <- "data/brca_xenvis/xenium"
data <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
data <- subset(data, subset = nCount_Xenium > 0)

VlnPlot(data, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

data <- SCTransform(data, assay = "Xenium")
data <- RunPCA(data, npcs = 30, features = rownames(data))
data <- RunUMAP(data, dims = 1:30)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.3)

p1 <- DimPlot(data, raster = FALSE)
p2 <- ImageDimPlot(data)
p1 + p2

save(data, file = "brca_xenium.rds")
