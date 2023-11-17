library(caret)
library(plyr)

#brca 10x
load("rds/brca_10x_major_all_annots.rds")
df <- data@meta.data %>% select(patho_annots, customclassif, rctd_annots, spotlight_annots)
load("rds/brca_ref_markers_major_normalized.rds")
from <- ref_markers$cluster %>% as.character() %>% unique()
to <- c(NA, "Fibrous Tissue", NA, "Immune Cells", "Immune Cells", NA, NA, "Immune Cells", "Invasive Carcinoma")
df$truth <- na_if(df$patho_annot, "Fat") %>% na_if("Necrosis")
df$sctype <- mapvalues(df$customclassif, from, to)
df$rctd <- mapvalues(df$rctd_annots, from, to)
df$spotlight <- mapvalues(df$spotlight_annots, from, to)
sctype_cm <- confusionMatrix(as.factor(df$sctype), as.factor(df$truth), mode = "everything")
rctd_cm <- confusionMatrix(as.factor(df$rctd), as.factor(df$truth), mode = "everything")
spotlight_cm <- confusionMatrix(as.factor(df$spotlight), as.factor(df$truth), mode = "everything")

acc <- c(Sctype = sctype_cm$overall[1],
         RCTD = rctd_cm$overall[1],
         SPOTlight = spotlight_cm$overall[1])
plot(acc)
axis(1, at = 1:3, c("Sctype", "RCTD", "SPOTlight"))


#brca xenvis visium
load("rds/brca_xenvis_visium.rds")
df <- data@meta.data %>% select(truth, customclassif, seurat_annots, rctd_annots, spotlight_annots)
load("rds/brca_ref2_markers.rds")

from <- ref_markers$cluster %>% as.character() %>% unique()
from
data@meta.data$truth %>% unique()
to <- c("stromal", NA, NA, "myoepithelial", "immune", "DCIS #1", "invasive", "invasive", "immune",
        "endothelial", NA, "DCIS #2", "immune", "stromal|immune", "immune|invasive", NA, NA, NA,
        "myoepithelial")

cm <- function(annots, from, to, truth) {
  mapped <- mapvalues(df[[annots]], from, to)
  mask <- mapply(grepl, pattern = mapped, x = truth)
  tr <- ifelse(mask, mapped, truth)
  tr <- ifelse(tr %in% mapped, tr, NA)
  res <- confusionMatrix(as.factor(mapped), as.factor(tr), mode = "everything")
  return(list(res, res$overall[1]))
}
plot_cm <- function(cm) {
  plt <- as.data.frame(cm$table)
  plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
  
  ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
    geom_tile() + geom_text(aes(label=Freq)) +
    scale_fill_gradient(low="white", high="#009194") +
    labs(x = "Reference",y = "Prediction") +
    scale_x_discrete() +
    scale_y_discrete()
}

sctype_cm <- cm("customclassif", from, to, df$truth)
seurat_cm <- cm("seurat_annots", from, to, df$truth)
rctd_cm <- cm("rctd_annots", from, to, df$truth)
spotlight_cm <- cm("spotlight_annots", from, to, df$truth)
sctype_cm[[2]]
rctd_cm[[2]]
seurat_cm[[2]]
spotlight_cm[[2]]
plot_cm(sctype_cm[[1]])
plot_cm(rctd_cm[[1]])
plot_cm(seurat_cm[[1]])
plot_cm(spotlight_cm[[1]])

accuracies = data.frame(method = c("sctype", "seurat", "rctd", "spotlight"),
                        accuracy = c(sctype_cm[[2]], seurat_cm[[2]],
                                     rctd_cm[[2]], spotlight_cm[[2]]))
ggplot(accuracies, aes(x=factor(method, levels = c("sctype", "rctd", "spotlight", "seurat")), y=accuracy)) + 
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  xlab("method")

#brca xenvis xenium
load("rds/brca_xenvis_xenium.rds")
df <- data@meta.data %>% select(truth, customclassif, rctd_annots, seurat_annots, spotlight_annots)
load("rds/brca_ref2_markers.rds")
from <- ref_markers$cluster %>% as.character() %>% unique()
from
data@meta.data$truth %>% unique()
to <- c("Stromal", "Macrophage 1", NA, "ACTA2+ Myoepithelial", "T Cells", "DCIS", "Invasive Tumor",
        "Invasive Tumor", "T Cells", "Endothelial", "Macrophage 2", "DCIS", "B Cells", "Stromal|T Cells",
        "T Cells|Invasive Tumor", NA, "Mast Cells", NA, "KRT15+ Myoepithelial")

cmx <- function(annots, from, to, truth) {
  mapped <- mapvalues(df[[annots]], from, to)
  mask <- mapply(grepl, pattern = truth, x = mapped)
  tr <- ifelse(mask, mapped, truth)
  tr <- ifelse(tr %in% mapped, tr, NA)
  res <- confusionMatrix(as.factor(mapped), as.factor(tr), mode = "everything")
  return(list(res, res$overall[1]))
}

sctype_cm <- cmx("customclassif", from, to, df$truth)
seurat_cm <- cmx("seurat_annots", from, to, df$truth)
rctd_cm <- cmx("rctd_annots", from, to, df$truth)
spotlight_cm <- cmx("spotlight_annots", from, to, df$truth)
sctype_cm[[2]]
rctd_cm[[2]]
seurat_cm[[2]]
spotlight_cm[[2]]
plot_cm(sctype_cm[[1]])
plot_cm(rctd_cm[[1]])
plot_cm(seurat_cm[[1]])
plot_cm(spotlight_cm[[1]])

ggplot(times, aes(x=method, y=mean)) + 
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  geom_errorbar( aes(x=method, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)

