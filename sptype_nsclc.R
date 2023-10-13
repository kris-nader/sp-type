library(Seurat)

#Annotate the data by integrating with single cell reference data
#read reference
#(ref from https://www.nature.com/articles/s41597-023-02074-6)
ref <- readRDS("nsclc_ref/RNA_rawcounts_matrix.rds")
metadata <- read.csv("nsclc_ref/metadata.csv")
ref <- CreateSeuratObject(ref)
ref@meta.data$celltype_major <- metadata$Cell_Cluster_level1
ref@meta.data$celltype_minor <- metadata$Cell_Cluster_level2
ref@meta.data$subtype <- metadata$Subtype

#set minor (or major) cell types as main identifier
ref <- SetIdent(ref, value = ref@meta.data$celltype_major)

ref <- subset(ref, subtype == "squamous cell carcinoma")
ref <- SCTransform(ref)
ref <- RunPCA(ref)
ref <- RunUMAP(ref, dims = 1:30)

#find markers of reference and save
ref_markers <- FindAllMarkers(reference)
save(ref_markers, file = "nsclc_scc_ref_markers_major_normalized.rds")

#you can also use second ref (nsclc_ref2.rds) from https://www.nature.com/articles/s41591-018-0096-5
load("nsclc_ref2.rds")

#load reference cluster markers
load("nsclc_ref_markers_major.rds")
load("nsclc.rds")

# I did not try this yet because couldn't normalize the reference
retrieve_markers(data, ref_markers, tissue="Lung")
data = run_plot_sctype(data = data, db_ = "temp/sctypeDB_20.xlsx", tissue = "Lung")
