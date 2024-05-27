library(SeuratData)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)


## reproduce fig 1B

options(timeout=600)
InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.5, verbose = FALSE)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2

source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");
slide.seq=run_sctype(slide.seq,known_tissue_type = "Hippo", slot = "SCT", name = "sptype_classification",custom_marker_file = "https://raw.githubusercontent.com/kris-nader/sp-type/master/ref_markers_hippo.xlsx")
my_colors <- list(
  "Astrocyte" = "#D6272A",
  "CA1 Principal cells" = "#1F78B4",
  "CA3 Principal cells" = "#A6CEE3",
  "Entorhinal cortex" = "#FB9A99",
  "Oligodendrocyte" = "#B2DF8A",
  "Ependymal" = "#33A02C",
  "Microglia" = "#FFFF99",
  "Interneuron" = "#FB820E",
  "Dentate Principal cells" = "#CAB2D6",
  "Endothelial tip"= "#B15928"
)

 SpatialDimPlot(slide.seq, stroke = 0,group.by="sptype_classification",cols=my_colors)


 ## reproduce fig 1C 
 ## anterior
 source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");
 load("/Users/naderkri/Desktop/sptype/brain/mouse_brain_anterior_10x.rds")
 brain_anterior=run_sctype(data,known_tissue_type = "Brain",slot="SCT",name="sptype_classification")
 my_brain=list(
   "Astrocytes" = "#D6272A",
   "Dopaminergic neurons" = "#B2DF8A",
   "GABAergic neurons" = "#34A02C",
   "Myelinating Schwann cells"= "#FB9A99",
   "Neuroblasts"="#FFFF99",
   "Glutamatergic neurons"= "#FF7F02",
   "Immature neurons"="#A6CEE3",
   "Mature neurons"="#1F78B4",
   "Radial glial cells"="#CAB2D6",
   "Oligodendrocytes"="#6A3D9A",
   "Endothelial cells"="#B15928"
 )
 
SpatialDimPlot(brain_anterior,group.by ="sptype_classification",cols=my_brain)
 
## posterior
load("/Users/naderkri/Desktop/sptype/brain/mouse_brain_posterior_10x.rds")
brain_posterior=run_sctype(data,known_tissue_type = "Brain",slot="SCT",name="sptype_classification")
SpatialDimPlot(brain_posterior,group.by ="sptype_classification",cols=my_brain)
SpatialDimPlot(brain_posterior,group.by ="sptype_classification",cols=my_brain)
 
## visium
load("/Users/naderkri/Desktop/sptype/visium_breast/brca_xenvis_visium.rds")
source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");

db_ <- "/Users/naderkri/Desktop/sptype/visium_breast/reference_datasets/ref_markers_breast_rm_mast_merged_frp_20_default_04.xlsx";
 
breast=run_sctype(data,known_tissue_type = "Breast",custom_marker_file = db_,slot="SCT",name="sptype_classification")
breast_cols <- list(
   "Adipocytes" = "#D6272A",
   "CD4+ T Cells" = "#1F78B4",
   "CD8+ T Cells" = "#1F78B4",
   "Dendritic cells" = "#A6CEE3",
   "DCs" = "#A6CEE3",
   "B cells" = "#8ED3C7",
   "B Cells" = "#8ED3C7",
   "Macrophages 1" = "#6A3D9A",
   "Macrophages 2" = "#6A3D9A",
   "DCIS 1" = "#B2DF8A",
   "DCIS 2" = "#34A02C",
   "Invasive Tumor" = "#006E2C",
   "Myoepi ACTA2+" = "#B15928",
   "Myoepi KRT15+"="#B15928",
   "Stromal"= "#FB9A99",
   "Stromal & T Cell Hybrid"="#FB9A99",
   "Endothelial"= "#FF7F02"
 )
 
SpatialDimPlot(breast,group.by ="sptype_classification",cols = breast_cols)
SpatialDimPlot(breast,group.by ="truth")

 ## xenium
load("/Users/naderkri/Desktop/sptype/visium_breast/brca_xenvis_xenium.rds")
source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");

db_ <- "/Users/naderkri/Desktop/sptype/visium_breast/reference_datasets/ref_markers_breast_rm_mast_merged_frp_20_default_04.xlsx";
breast=run_sctype(data,known_tissue_type = "Breast",custom_marker_file = db_,slot="SCT",name="sptype_classification")
ImageDimPlot(breast,group.by ="sptype_classification",cols = breast_cols)


## reproduce fig 1E
load("/Users/naderkri/Desktop/sptype/visium_breast/brca_xenvis_visium.rds")

df <- breast@meta.data %>% select(truth, sptype_classification, seurat_annots, rctd_annots,spotlight_annots)

from=unique(c(df$sptype_classification,df$seurat_annots,df$rctd_annots,df$spotlight_annots))
to=c("stromal/endothelial","immune","myoepithelial","immune","invasive",
     "stromal/endothelial","DCIS #1","DCIS #2","stromal/endothelial","immune",
     "myoepithelial","immune","invasive","stromal/endothelial","stromal/endothelial",
     "invasive","immune","immune","immune"
     )


from_truth=df$truth %>% unique()
to_truth=c("stromal/endothelial","mixed","stromal/endothelial","stromal/endothelial","invasive","invasive","stromal/endothelial","DCIS #1","immune","DCIS #2", "myoepithelial")

df$truth=mapvalues(df$truth,from=from_truth,to=to_truth)

df$sptype_classification=mapvalues(df$sptype_classification,from = from,to=to)
df$seurat_annots=mapvalues(df$seurat_annots,from = from,to=to)
df$rctd_annots=mapvalues(df$rctd_annots,from = from,to=to)
df$spotlight_annots=mapvalues(df$spotlight_annots,from = from,to=to)


df_sub=df %>% filter(truth!="mixed")


sctype_cm <- confusionMatrix(as.factor(df_sub$sptype_classification), factor(df_sub$truth,levels=unique(df_sub$truth)), mode = "everything")
rctd_cm <- confusionMatrix(as.factor(df_sub$rctd_annots), factor(df_sub$truth,levels=unique(df_sub$truth)), mode = "everything")
spotlight_cm <- confusionMatrix(as.factor(df_sub$spotlight_annots), factor(df_sub$truth,levels=unique(df_sub$truth)), mode = "everything")
seurat_cm <- confusionMatrix(as.factor(df_sub$seurat_annots), factor(df_sub$truth,levels=unique(df_sub$truth)), mode = "everything")

df=data.frame(ScType=sctype_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"F1"],
              Seurat=seurat_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"F1"],
              RCTD=rctd_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"F1"],
              SPOTlight=spotlight_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"F1"])

df1=data.frame(ScType=sctype_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Precision"],
              Seurat=seurat_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Precision"],
              RCTD=rctd_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Precision"],
              SPOTlight=spotlight_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Precision"])


df2=data.frame(ScType=sctype_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Recall"],
              Seurat=seurat_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Recall"],
              RCTD=rctd_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Recall"],
              SPOTlight=spotlight_cm$byClass[paste0("Class: ",unique(df_sub$truth)),"Recall"])





rownames(df)=c("stromal/endothelial","invasive","DCIS 1","immune","DCIS 2","myoepithelial")
rownames(df1)=c("stromal/endothelial","invasive","DCIS 1","immune","DCIS 2","myoepithelial")
rownames(df2)=c("stromal/endothelial","invasive","DCIS 1","immune","DCIS 2","myoepithelial")



df$cell_types=rownames(df)
df1$cell_types=rownames(df1)
df2$cell_types=rownames(df2)


# load tidyverse 
library(tidyverse) 

# Convert to long format
df_long <- df %>%
  pivot_longer(cols = -cell_types, names_to = "tool", values_to = "F1")


df_long1 <- df1 %>%
  pivot_longer(cols = -cell_types, names_to = "tool", values_to = "Precision")


df_long2 <- df2 %>%
  pivot_longer(cols = -cell_types, names_to = "tool", values_to = "Recall")

# Define shapes for each cell type
shapes <- c("stromal/endothelial" = 17, "invasive" = 18,"DCIS 1" = 15,"DCIS 2"=12,"myoepithelial"=13,"immune"=10)

# Create the plot
pdf("/Users/naderkri/Desktop/sptype/benchmark/visium_F1_cell_types.pdf",height=3,width = 3)
p <- ggplot(df_long, aes(x = factor(tool, levels=c("ScType","Seurat","RCTD","SPOTlight")), y = F1, group = cell_types, color = cell_types, shape = cell_types)) +
  geom_line() +
  geom_point(size = 3) +
  scale_shape_manual(values = shapes) +
  theme_classic() +
  labs(x = "Tool", y = "F1", color = "Cell Type", shape = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("ScType"="#E21A1B", "Seurat"="#FDBE6E", "RCTD"="#2078B4", "SPOTlight"="#A6CEE3"))+
  ylim(0,1)+NoLegend()
p
dev.off()

pdf("/Users/naderkri/Desktop/sptype/benchmark/visium_precision_cell_types.pdf",height=3,width = 3)
p1 <- ggplot(df_long1, aes(x = factor(tool, levels=c("ScType","Seurat","RCTD","SPOTlight")), y = Precision, group = cell_types, color = cell_types, shape = cell_types)) +
  geom_line() +
  geom_point(size = 3) +
  scale_shape_manual(values = shapes) +
  theme_classic() +
  labs(x = "Tool", y = "Precision", color = "Cell Type", shape = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("ScType"="#E21A1B", "Seurat"="#FDBE6E", "RCTD"="#2078B4", "SPOTlight"="#A6CEE3"))+
  ylim(0,1)+NoLegend()
p1
dev.off()

pdf("/Users/naderkri/Desktop/sptype/benchmark/visium_recall_cell_types.pdf",height=3,width = 3)
p2 <- ggplot(df_long2, aes(x = factor(tool, levels=c("ScType","Seurat","RCTD","SPOTlight")), y = Recall, group = cell_types, color = cell_types, shape = cell_types)) +
  geom_line() +
  geom_point(size = 3) +
  scale_shape_manual(values = shapes) +
  theme_classic() +
  labs(x = "Tool", y = "Recall", color = "Cell Type", shape = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("ScType"="#E21A1B", "Seurat"="#FDBE6E", "RCTD"="#2078B4", "SPOTlight"="#A6CEE3"))+
  ylim(0,1)+NoLegend()
p2
dev.off()

p+p1+p2 

# Plot
df_long %>%
  ggplot( aes(x=factor(tool, levels=c("ScType","Seurat","RCTD","SPOTlight")), y=F1, fill=tool)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black")+
  scale_fill_manual(values=c("ScType"="#E21A1B", "Seurat"="#FDBE6E", "RCTD"="#2078B4", "SPOTlight"="#A6CEE3"))+
  theme_classic()

df_long1 %>%
  ggplot( aes(x=factor(tool, levels=c("ScType","Seurat","RCTD","SPOTlight")), y=Precision, fill=tool)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black")+
  scale_fill_manual(values=c("ScType"="#E21A1B", "Seurat"="#FDBE6E", "RCTD"="#2078B4", "SPOTlight"="#A6CEE3"))+
  theme_classic()

df_long2 %>%
  ggplot( aes(x=factor(tool, levels=c("ScType","Seurat","RCTD","SPOTlight")), y=Recall, fill=tool)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black")+
  scale_fill_manual(values=c("ScType"="#E21A1B", "Seurat"="#FDBE6E", "RCTD"="#2078B4", "SPOTlight"="#A6CEE3"))+
  theme_classic()



 ## supplement fig 
 load("/Users/naderkri/Desktop/sptype/visium_breast/brca_xenvis_visium.rds")
 #load("/Users/naderkri/Desktop/sptype/visium_breast/brca_xenvis_xenium.rds")
 my_visium <- list(
   "DCIS #1" = "#A6CEE3",
   "DCIS #2" = "#1F78B4",
   "adipocytes" = "#E31A1C",
   "immune" = "#B2DF8A",
   "invasive" = "#FB9A99",
   "mixed"="#FDBF6F",
   "mixed/invasive"="#CAB2D6",
   "myoepithelial/stromal/immune"="#33A02C",
   "stromal"="#FFFF99",
   "stromal/endothelial"="#6A3D9A",
   "stromal/endothelial/immune"="#B15928"
 )
 my_xenium <- list(
   "DCIS" = "#A6CEE3",
   "T Cells" = "#E5F5E0",
   "Macrophage 1" = "#C7E9C0",
   "Macrophage 2" = "#A1D99B",
   "NK Cells"="#74C476",
   "Plasmablast"="#41AB5D",
   "B Cells"="#238B45",
   "Mast Cells"="#006D2C",
   "Plasmacytoid Dendritic"="#00441B",
   "ACTA2+ Myoepithelial"="#B15928",
   "KRT15+ Myoepithelial"="#FB9A99",
   "Invasive Tumor" = "#1F78B4",
   "Stromal"="#FFFF99",
   "Endothelial"="#E41A1C",
   "Undefined"="lightgray",
   "<NA>"="lightgray"
 )
 
 ggplot(ggData, aes(sptype, value, fill = truth)) +
   geom_col() + xlab("sptype") + ylab("Proportion of Cells (%)") + coord_flip()+theme_classic()+
   scale_fill_manual(values = my_visium)
 
 
 
