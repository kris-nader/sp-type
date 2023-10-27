library(Seurat)
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

ref <- SCTransform(ref)
ref <- RunPCA(ref)
ref <- RunUMAP(ref, dims = 1:30)

#set minor (or major) cell types as main identifier
ref <- SetIdent(ref, value = ref@meta.data$celltype_minor)

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

source("sptype_KMN.R")
load("rds/brca_xenvis_xenium.rds")
load("rds/brca_ref2_markers.rds")

marker_per_cluster <- NULL
tissue <- "Breast"
markers <- retrieve_markers(ref_markers,
                 tissue = tissue,
                 marker_per_cluster = marker_per_cluster,
                 min_pct_diff = 0.40)
data <- run_plot_sctype(data = data,
                db_ = paste0("temp/sctypeDB_", marker_per_cluster, ".xlsx"),
                tissue = tissue,
                marker_per_cluster = marker_per_cluster,
                output_folder = "figures/brca_xenvis/xenium",
                output_name = "sptype_40",
                pt.size.factor = 2.4,
                saveRDS = FALSE,
                st_method = "xenium",
                )
save(data, file = "rds/brca_xenvis_xenium.rds")
#check if retrieved markers have similar pct.1 and pct.2
a <- ref_markers[abs(ref_markers$pct.1 - ref_markers$pct.2) < 0.1,]
m <- c(markers$geneSymbolmore1 %>% strsplit(",") %>% unlist(),
       markers$geneSymbolmore2 %>% strsplit(",") %>% unlist()) %>% unique()
genes_to_check <- a[a$gene %in% m,]
bad_genes <- genes_to_check$gene %>% unique

#0.2 564/696
#0.1 358/696
#0.05 223/696

