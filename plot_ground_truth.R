library(openxlsx)

load("rds/brca_xenvis_visium.rds")
annots <- read.xlsx("data/cancer/brca_xenvis/Cell_Barcode_Type_Matrices.xlsx", "Visium")
barcodes <- data.frame(Barcode = data@assays$Spatial@counts@Dimnames[[2]])
annots <- merge(annots, barcodes, all.y = TRUE, by = "Barcode")
data@meta.data$truth <- annots$Annotation
SpatialDimPlot(data, group.by = "truth", pt.size.factor = 2.4)
save(data, file = "rds/brca_xenvis_visium.rds")

load("rds/brca_xenvis_xenium.rds")
annots <- read.xlsx("data/cancer/brca_xenvis/Cell_Barcode_Type_Matrices.xlsx", "Xenium R1 Fig 3 (unsupervised)")
barcodes <- data.frame(Barcode = data@assays$Xenium@counts@Dimnames[[2]])
annots <- merge(annots, barcodes, all.y = TRUE, by = "Barcode")
data@meta.data$truth <- annots$ident
ImageDimPlot(data, group.by = "truth")
save(data, file = "rds/brca_xenvis_xenium.rds")
