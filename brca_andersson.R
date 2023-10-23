library(Seurat)
library(data.table)
library(ggplot2)
library(plotly)
library(zeallot)
library(openxlsx)
library(STutility)

meta_data <- read.xlsx("data/brca_andersson/clinical_data/10_HER2+_info.xlsx")
rownames(meta_data) <- meta_data$Sample
samples <- list.files(pattern = ".tsv", path = "data/brca_andersson/ST-cnts/", full.names = T)
names(samples) <- substr(do.call(rbind, strsplit(samples, split = "/"))[, 4], start = 1, stop = 2)
imgs <- list.files(path = "data/brca_andersson/ST-imgs/", recursive = T, full.names = T, pattern = ".jpg")
names(imgs) <- do.call(rbind, strsplit(imgs, split = "/"))[, 5]
ids <- names(samples)
infoTable <- data.frame(samples, imgs = imgs[ids], ids, patient_id = substr(x = ids, start = 1, stop = 1), stringsAsFactors = FALSE)
infoTable <- cbind(infoTable, meta_data[infoTable$patient_id, ])
infoTable[, 8:ncol(infoTable)]

infoTable$spotfiles <- list.files(path = "data/brca_andersson/ST-spotfiles", full.names = T)[1:36]
head(infoTable)

seu.list <- lapply(unique(infoTable$ids), function(s) {
  InputFromTable(infotable = subset(infoTable, ids == s), 
                 min.gene.spots = 20,
                 min.spot.feature.count = 300,
                 platform = "1k")
}) 

seu.list <- lapply(seu.list, function(seu) {
  seu %>% LoadImages(verbose = T, time.resolve = F)
})

#choose one section
data <- seu.list[[1]]
FeatureOverlay(data, features = "nFeature_RNA")

#filter spots with high mitochondrial genes (min feature count filtering done above)
mt.genes <- grep(pattern = "^mt-", x = rownames(data), value = TRUE)
data$percent.mito <- (Matrix::colSums(data@assays$RNA@counts[mt.genes, ])/Matrix::colSums(data@assays$RNA@counts))*100
FeatureOverlay(data, features = "percent.mito")
data <- SubsetSTData(data, expression = percent.mito < 30)
cat("Spots removed: ", ncol(data) - ncol(data.subset), "\n")

data <- SCTransform(data)
data <- RunPCA(data, assay = "SCT", verbose = TRUE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.8) #0.8 is the default resolution
data <- RunUMAP(data, reduction = "pca", dim = 1:30)
save(data, file = "brca_andersson_A1.rds")
