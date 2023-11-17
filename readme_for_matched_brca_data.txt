sctype:
run sptype_brca.R lines 42-62
data: homes/knader/misra/ST_fimm/rds/brca_xenvis_visium.rds
markers: homes/knader/misra/ST_fimm/rds/brca_ref2_markers.rds

seurat:
run brca_seurat_annot.R lines 1-9 to scttransform the data, then after 18
data: homes/knader/misra/ST_fimm/rds/brca_xenvis_visium.rds
ref: homes/knader/misra/ST_fimm/rds/brca_ref2_normalized.rds

rctd:
run brcac_rctd.R lines 15-42
data: homes/knader/misra/ST_fimm/rds/brca_xenvis_visium.rds
ref: homes/knader/misra/ST_fimm/data/cancer/brca_xenvis/sc/ftp
metadata: homes/knader/misra/ST_fimm/data/cancer/brca_xenvis/Cell_Barcode_Type_Matrices.xlsx
run below 42 to add annotations to the seurat object and for some plots

spotlight:
run brca_spotlight.R lines 1-67
data: homes/knader/misra/ST_fimm/rds/brca_xenvis_visium.rds
ref: homes/knader/misra/ST_fimm/rds/brca_ref2_normalized.rds
run below 67 to add annotations to the seurat object and for some plots

