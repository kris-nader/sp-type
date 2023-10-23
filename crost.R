library(Seurat)

move_file <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}


tissue <- "kidney"
dataset_num <- "000029"
sample_num <- "000134"
data_dir <- paste0("data/healthy/", tissue, "_crost_VISDP", dataset_num, "/VISDS", sample_num)
dir.create(data_dir, recursive = TRUE)
filename <- paste0("download_page_VISDS", sample_num, "_filtered_feature_bc_matrix.h5")
move_file(paste0("C:/Users/misratasci/Downloads/", filename), paste0(data_dir, "/", filename))
data <- Load10X_Spatial(data.dir = data_dir, filename = filename)

#quality control
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
plot1 <- VlnPlot(data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                 pt.size = 0.1)
plot2 <- SpatialFeaturePlot(data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                            pt.size.factor = 2, alpha = c(0.3, 1))

#save quality control plots
output_folder <- paste0("figures/crost/", tissue, "/VISDP", dataset_num, "/VISDS", sample_num)
dir.create(output_folder, recursive = TRUE)
output_name <- paste0(tissue, "_VISDS", sample_num, "_violin")
png(file = paste0("./", output_folder, "/", output_name, ".png"),
    width=1000, height=550)
print(plot1)
dev.off()
output_name <- paste0(tissue, "_VISDS", sample_num, "_quality_control")
png(file = paste0("./", output_folder, "/", output_name, ".png"),
    width=1000, height=550)
print(plot2)
dev.off()
plot1
#filter
data <- data[, data$nFeature_Spatial > 300 & data$percent.mt < 40]

#normalization
data <- SCTransform(data, assay="Spatial")

#dim reduction and clustering
data <- RunPCA(data, assay = "SCT", verbose = TRUE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.8) #0.8 is the default resolution
data <- RunUMAP(data, reduction = "pca", dim = 1:30)

#sptype
source("sptype_KMN.R")
library(stringr)
ptsize = 2
data <- run_plot_sctype(data = data,
                        #db_ = "temp/sctypeDB_20.xlsx",
                        tissue = str_to_title(tissue),
                        output_folder = output_folder,
                        output_name = paste0(tissue, "_crost_VISDS", sample_num),
                        saveRDS = TRUE,
                        pt.size.factor = ptsize)

#add crost's celltypes to data
celltype_deco <- read.csv(paste0(data_dir, "/download_page_VISDS", sample_num, "_celltype_deco.csv"),
                          row.names = 1)
#for each barcode get the cell type with highest weight
celltype_deco$celltype <- colnames(celltype_deco)[max.col(celltype_deco, ties.method = "random")]
barcodes <- data.frame(X = data@assays$Spatial@counts@Dimnames[[2]])
celltype_deco$X = row.names(celltype_deco)
df <- merge(celltype_deco, barcodes, all.y = TRUE, by = "X")
data@meta.data$crost_celltype_deco <- df$celltype

#plot crost's deco
p1 <- DimPlot(data, reduction = "umap", label = FALSE, repel = TRUE, group.by = 'crost_celltype_deco') 
p2 <- SpatialDimPlot(data, group.by = 'crost_celltype_deco', pt.size.factor = ptsize, alpha = c(1, 1), label = FALSE)
output_name <- paste0(tissue, "_VISDS", sample_num, "_crost_celltype_deco")
png(file = paste0("./", output_folder, "/", output_name, ".png"),
    width=1000, height=550)
print(p1 + p2)
dev.off()

#add additional metadata (for liver dataset 66)
#the number of barcodes in metadata and in data don't match
# metadata <- read.csv("data/healthy/liver_crost_VISDP000066/paper_annot.csv")
# s350 <- metadata[metadata$sample == "JBO14",]
# s351 <- metadata[metadata$sample == "JBO15",]
# s352 <- metadata[metadata$sample == "JBO18",]
# s353 <- metadata[metadata$sample == "JBO19",]
# s354 <- metadata[metadata$sample == "JBO22",]
# data <- readRDS("liver_crost_VISDS000350.RDS")
# barcodes <- data.frame(X = data@assays$Spatial@counts@Dimnames[[2]])
# s350$spot = substr(s350$spot,1,nchar(s350$spot)-2)
# s350$X <- s350$spot
# s350 <- subset(s350, select = c(X, zonationGroup))
# df <- merge(s350, barcodes, all.y = TRUE, by = "X")
# data@meta.data$paper_annot <- df$zonationGroup
# SpatialDimPlot(data, group.by = "paper_annot")
