library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#load data
data_dir <- "data/brca_xenvis/visium"
filename <- "filtered_feature_bc_matrix.h5"
data <- Load10X_Spatial(data.dir = data_dir, filename = filename)

#quality control
#nFeature is the number of genes detected in each cell
#nCount is the total number of molecules detected within a cell
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
plot1 <- VlnPlot(data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                 pt.size = 0.1)
plot2 <- SpatialFeaturePlot(data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                            pt.size.factor = 2.4, alpha = c(0.3, 1))

#filter
data <- data[, data$nFeature_Spatial > 500 & data$percent.mt < 20]

#see top expressed genes
C <- data@assays$Spatial@counts
C@x <- C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

#normalization
data <- SCTransform(data, assay="Spatial")
#without using SCT:
#data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000) %>% ScaleData() %>% FindVariableFeatures()

#spatial plot for one gene
#SpatialFeaturePlot(data, features = c("BRCA1"), pt.size.factor = 2.4, alpha = c(0.1, 1))

#dim reduction and clustering
data <- RunPCA(data, assay = "SCT", verbose = TRUE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.8) #0.8 is the default resolution
data <- RunUMAP(data, reduction = "pca", dim = 1:30)

#save data as RData
save(data, file = "brca_xenvis_visium.RData")

#umap plot
p1 <- DimPlot(data, reduction = "umap", label = FALSE)
#spatial umap plot
p2 <- SpatialDimPlot(data, label = FALSE, label.size = 3, pt.size.factor = 2.4, alpha = c(1, 1))
p1 + p2
#plot clusters seperately
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