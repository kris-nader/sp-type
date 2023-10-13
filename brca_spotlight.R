library(ggplot2)
library(SPOTlight)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)

load("brca.rds")
load("brca_ref.rds")
#ref <- subset(ref, subtype == "HER2+")

sce <- as.SingleCellExperiment(ref)
sce <- logNormCounts(sce)

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))

dec <- modelGeneVar(sce, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)

colLabels(sce) <- colData(sce)$celltype_minor

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$celltype_minor)
# downsample to at most 20 per identity & subset
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

res <- SPOTlight(
  x = sce,
  y = data,
  groups = as.character(sce$celltype_minor),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

save(res, file = "brca_minor_spotlight_res_all_ref.rds")

head(mat <- res$mat)[, seq_len(3)]

mod <- res$NMF

plotTopicProfiles(
  x = mod,
  y = sce$celltype_minor,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

plotTopicProfiles(
  x = mod,
  y = sce$celltype_minor,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)

library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)

plotCorrelationMatrix(mat)
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")

#scatterpie
ct <- colnames(mat)
mat[mat < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

plotSpatialScatterpie(
  x = data,
  y = mat,
  cell_types = colnames(mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))

plotSpatialScatterpie(
  x = data,
  y = mat,
  cell_types = colnames(mat),
  img = TRUE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))
