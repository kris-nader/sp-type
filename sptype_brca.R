library(Seurat)
source("sptype_KMN.R")
#Prepare reference

#read reference sc file
#ref from https://www.nature.com/articles/s41588-021-00911-1
ref <- Read10X("brca_ref", gene.column = 1)
ref <- CreateSeuratObject(ref)

#read metadata
metadata <- read.csv("brca_ref/metadata.csv")

#add metadata to reference
ref@meta.data$X <- metadata$X
ref@meta.data$percent.mito <- metadata$percent.mito
ref@meta.data$subtype <- metadata$subtype
ref@meta.data$celltype_subset <- metadata$celltype_subset
ref@meta.data$celltype_minor <- metadata$celltype_minor
ref@meta.data$celltype_major <- metadata$celltype_major

#set minor (or major) cell types as main identifier
ref <- SetIdent(ref, value = ref@meta.data$celltype_minor)

#added this part recently (10.11.23), data from before was made with non-normalized reference
ref <- subset(ref, subtype == "HER2+") # cannot normalize all of the data because weak computer :(
ref <- SCTransform(ref)
ref <- RunPCA(ref)
ref <- RunUMAP(ref, dims = 1:30)

#find markers for cell types
ref_markers <- FindAllMarkers(ref)
#save before filtering
save(ref_markers, file = "brca_her2_ref_markers_minor_normalized.rds")

#PVL: perivascular type cells

#add pathologist annotations to data
patho <- read.csv("brca_patho_annot.csv")
spots <- data@assays$SCT@counts@Dimnames[[2]]
spots <- data.frame(Barcode = spots)
patho <- join(patho, spots, by = "Barcode", type = "right")
data@meta.data$patho_annot <- patho$Pathologist.Annotation
save(data, file = "brca.rds")

load("brca.rds")
load("brca_her2_ref_markers_major_normalized.rds")

retrieve_markers(data, ref_markers, tissue="Breast")
data = run_plot_sctype(data = data, db_ = "temp/sctypeDB_20.xlsx", tissue = "Breast")




