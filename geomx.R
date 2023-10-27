library(GeomxTools)
library(Seurat)

dcc_files <- dir("data/nsclc_geomx/GSE200563_RAW", full.names = TRUE)
pkc_files <- dir("data/nsclc_geomx/pkc", full.names = TRUE)

data <- readNanoStringGeoMxSet(dcc_files, pkc_files)
