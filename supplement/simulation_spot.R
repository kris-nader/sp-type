#### process single cell data
lapply(c("dplyr","Seurat","HGNChelper","packcircles","plotrix","scales","plyr","reshape2",
         "matrixStats","scatterpie","cowplot","ggpubr","ggplot2","caret"), library, character.only = T)
## remove the mast cells, Stromal & T cells , T cell & Tumor hybrid 
## merge the Tumor with DCIS
'%ni%'=Negate('%in%')

set.seed(15)
reference=readRDS("/Users/naderkri/Desktop/sptype/visium_breast/scFFPE_done.RDS")
cells_=reference@meta.data[which(!is.na(reference@meta.data$annotations)),] %>% rownames
reference=subset(reference,cells = cells_)
## remove the stromal and t cell and tumor hybrids
cells_also=rownames(reference@meta.data[which(reference@meta.data$annotations %ni% 
                                                c("Mast Cells","Stromal & T Cell Hybrid","T Cell & Tumor Hybrid","Perivascular-Like",
                                                  "IRF7+ DCs","LAMP3+ DCs") ),])
reference=subset(reference,cells = cells_also)
table(reference@meta.data$annotations)

from_reference=as.character(reference$annotations %>% unique)
to_reference=c("Stromal","Macrophages","Myoepithelial",
               "T Cells","DCIS","DCIS","DCIS",
               "T Cells","Endothelial","Macrophages","DCIS",
               "B Cells","Myoepithelial")
reference$annotations_new=mapvalues(reference$annotations,from_reference,to_reference)

DimPlot(reference,group.by="annotations_new")
reference[["percent.mt"]] <- PercentageFeatureSet(reference, pattern = "^MT-")
Idents(reference)=1
VlnPlot(reference, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

reference <- subset(reference, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 6)

saveRDS(reference,file="/Users/naderkri/Desktop/sptype/benchmark/scffpe_filtered.RDS")

df <- data.frame(
  cell_name = rownames(reference@meta.data),
  annotation = reference@meta.data$annotations_new
)

# Set the seed for reproducibility
set.seed(123)
# Create a partition index
k <- createFolds(df$annotation,k = 2,list = F)
df$index=k

cell_names_training <- df %>%filter(index == 1) %>% pull(cell_name)
cell_names_testing <- df %>%filter(index == 2) %>% pull(cell_name)

training_seurat=subset(reference,cells = cell_names_training )
testing_seurat=subset(reference,cells = cell_names_testing )

saveRDS(training_seurat,"./training_seurat_1.RDS")
saveRDS(testing_seurat,"./testing_seurat_1.RDS")


reference=readRDS("./training_seurat_1.RDS")
DimPlot(reference, group.by="annotations_new")

reference <- SCTransform(reference, assay = "RNA", verbose = FALSE)
reference <- RunPCA(reference, assay = "SCT", verbose = FALSE) 
ElbowPlot(reference)
reference <- FindNeighbors(reference, reduction = "pca", dims = 1:50)
reference <- FindClusters(reference, verbose = FALSE,resolution = 0.8)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:50)

DimPlot(reference, group.by="annotations_new")

DimPlot(reference, group.by="seurat_clusters")
#cells_remove=rownames(reference@meta.data[which(reference@meta.data$seurat_clusters==14),])
cells_remove=CellSelector(DimPlot(reference, group.by="annotations_new"))
cells_keep=setdiff(Cells(reference),cells_remove)
reference=subset(reference,cells=cells_keep)

reference <- SCTransform(reference, assay = "RNA", verbose = FALSE)
reference <- RunPCA(reference, assay = "SCT", verbose = FALSE) 
ElbowPlot(reference)
reference <- FindNeighbors(reference, reduction = "pca", dims = 1:50)
reference <- FindClusters(reference, verbose = FALSE,resolution = 0.8)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:50)


DimPlot(reference, group.by="annotations_new")
saveRDS(reference,"./training_seurat_preprocessed_1.RDS")

testing=readRDS("./testing_seurat_1.RDS")
DimPlot(testing, group.by="annotations_new")

testing <- SCTransform(testing, assay = "RNA", verbose = FALSE)
testing <- RunPCA(testing, assay = "SCT", verbose = FALSE) 
ElbowPlot(testing)
testing <- FindNeighbors(testing, reduction = "pca", dims = 1:50)
testing <- FindClusters(testing, verbose = FALSE,resolution = 0.8)
testing <- RunUMAP(testing, reduction = "pca", dims = 1:50)
DimPlot(testing, group.by="annotations_new")

cells_remove=CellSelector(DimPlot(testing, group.by="annotations_new"))
cells_keep=setdiff(Cells(testing),cells_remove)
testing=subset(testing,cells=cells_keep)
DimPlot(testing, group.by="annotations_new")

testing <- SCTransform(testing, assay = "RNA", verbose = FALSE)
testing <- RunPCA(testing, assay = "SCT", verbose = FALSE) 
ElbowPlot(testing)
testing <- FindNeighbors(testing, reduction = "pca", dims = 1:50)
testing <- FindClusters(testing, verbose = FALSE,resolution = 0.8)
testing <- RunUMAP(testing, reduction = "pca", dims = 1:50)
DimPlot(testing, group.by="annotations_new")


saveRDS(testing,"./testing_seurat_preprocessed_1.RDS")


#####################################

reference=readRDS("./training_seurat_preprocessed_1.RDS")


## simulate the cells, cell types and sizes that will be used 
from=c(10/2,8/2,12/2,15/2,20/2)
number_single_cells=10000
#number_single_cells=1000
number_spots=6000
radius_spot=10
to=c("T Cells","B Cells","Stromal","DCIS","Macrophages")
color=c("#1f78b4","#33a02c","#ffff99","#fb9a99","#6a3d9a")

set.seed(15)

## here we imitate the cell type proportions from the original reference dataset 
areas <- sample(from, number_single_cells, rep = TRUE, prob = (unname(table(reference$annotations_new)[to])/sum(unname(table(reference$annotations_new)[to])))*100)

mapped=data.frame(area=areas,map=mapvalues(areas,from,to))
mapped$cell=""
for( i in 1:nrow(mapped) ){
  sampled.cell_C1 <- sample(x = reference@meta.data %>% filter(annotations_new==mapped[i,"map"]) %>% rownames(),size =1, replace = F)
  mapped$cell[i]=sampled.cell_C1 
}
packing <- data.frame(circleProgressiveLayout(areas,sizetype = "radius"))
packing$cell_name=mapped$cell
packing$cell_type=mapped$map
packing$color=mapvalues(packing$cell_type,to,color)
mapped=unique(mapped)
packing$both=paste0(packing$cell_name,"_",1:number_single_cells,"_",packing$cell_type)

areas_larger <- c(rep(radius_spot,number_spots))
packing_larger <- data.frame(circleProgressiveLayout(x=areas_larger,sizetype = "radius"))
packing_larger$spot_name=paste0("spot_",1:number_spots)


## Plot the circles-- for visualization
plot(seq(min(packing$x),max(packing$x),length=20),seq(min(packing$y),max(packing$y),length=20),type="n",xlab="",ylab="",main=paste0("Spatial Transcriptomics Simulation ",nrow(packing)," cells and ",nrow(packing_larger)," spots"))
for( z in 1:nrow(packing)){
  draw.circle(packing$x[z],packing$y[z],packing$radius[z],col =packing$color[z] )
}
for( y in 1:nrow(packing_larger)){
  draw.circle(x=packing_larger$x[y],packing_larger$y[y],packing_larger$radius[y],col=alpha("gray",0.15))
}


##################################
## function to calculate the distance between 2 centers of 2 circles
distance <- function(x1, y1, x2, y2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}


## function to calculate intersection area between two circles
circle_intersection <- function(x1, y1, r1, x2, y2, r2) {
  d <- distance(x1, y1, x2, y2)
  # Check if circles are completely separate
  if (d > r1 + r2) {
    return(0)
  }
  # Check if one circle is inside the other
  if (d + min(r1,r2) <= max(r1,r2)) {
    r <- min(r1, r2)
    return(pi * r^2)
  }
  # Calculate intersection area
  a1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  a2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  a3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
  return(a1 + a2 - a3)
}


spots_cells <- matrix(0, nrow = nrow(packing_larger), ncol = nrow(packing))
colnames(spots_cells)=packing$both
rownames(spots_cells)=packing_larger$spot_name

for (j in 1:nrow(packing_larger)){
  larger_circle <- data.frame(x = packing_larger$x[j], y = packing_larger$y[j], radius = packing_larger$radius[j],spot=packing_larger$spot_name[j])
  find_smaller_cells_around= packing %>% filter(abs(x) <= 5*(abs(larger_circle$x)+larger_circle$radius) &
                                                  abs(y) <= 5*(abs(larger_circle$y)+larger_circle$radius))
  #draw.circle(x=larger_circle$x,larger_circle$y,larger_circle$radius,col=alpha("blue",0.25))
  
  for (i in 1:nrow(find_smaller_cells_around)){
    smaller_circle <- data.frame(x = find_smaller_cells_around$x[i], y = find_smaller_cells_around$y[i], radius = find_smaller_cells_around$radius[i],cell=find_smaller_cells_around$cell_name[i])
    #draw.circle(x=smaller_circle$x,smaller_circle$y,smaller_circle$radius,col=alpha("red",0.35))
    smaller_circle_area <- pi * smaller_circle$radius^2
    intersection_area <- tryCatch({
      circle_intersection(x1 = larger_circle$x, y1 = larger_circle$y, r1 = larger_circle$radius,
                          x2 = smaller_circle$x, y2 = smaller_circle$y, r2 = smaller_circle$radius)
    }, warning = function(w) {
      return(0)
    })
    # Calculate the percentage of smaller circle area enclosed in larger circle
    percentage_enclosed <- (intersection_area / smaller_circle_area) * 100
    fix_index=as.integer(strsplit(find_smaller_cells_around[i,"both"],split = "_")[[1]][2])
    spots_cells[j, fix_index] <- ifelse((percentage_enclosed < 0| is.na(percentage_enclosed)), 0, percentage_enclosed)
  }
}



## remove spots with no cells and cells not covered by any spots
rm_these_cells=names(which(colSums(spots_cells)==0))
rm_these_spots=names(which(rowSums(spots_cells)==0))
spots_cells_rm=spots_cells
if(!identical(rm_these_cells,character(0))){
  spots_cells_rm <- spots_cells[, -which(colnames(spots_cells) %in% rm_these_cells)]
}
if(!identical(rm_these_spots,character(0))){
  spots_cells_rm <- spots_cells_rm[-which(rownames(spots_cells_rm) %in% rm_these_spots), ]
}



sparse_k=spots_cells_rm
split_data=colsplit(colnames(sparse_k), "\\_", names=c("A", "B"))

colnames(sparse_k)=mapvalues(split_data$A,from=mapped$cell,to=mapped$map)
sparse_mc_type=t(rowsum(t(sparse_k), group = colnames(sparse_k), na.rm = T))
sparse_mc_type_1=(sparse_mc_type/rowSums(sparse_k))*100

# filter spots where the max is less than 75-- in this way we make sure we have a ground truth
sparse_l=spots_cells_rm[names(which(apply(sparse_mc_type_1, 1, max)>=75)),]


remove_spots_clean_=names(which(is.na(rowSums(sparse_mc_type_1))))
if(!identical(remove_spots_clean_,character(0))){
  sparse_l=sparse_l[-which(rownames(sparse_l) %in% remove_spots_clean_),]
}

## now for each spot, aggregate counts from the reference dataset
merged3=matrix(0,nrow=dim(reference@assays$RNA$counts),ncol=0)
k=rownames(sparse_l)
for ( i in 1:nrow(sparse_l)){
  cells_index=colsplit(names(which(sparse_l[i, ] > 0)), "\\_", names=c("A", "B"))
  if(dim(cells_index)[1]>0){
    subset=reference[,c(cells_index$A)]
    subset$sctype_classification=k[i]
    if(length(Cells(subset)) > 1) {
      w=unname(sparse_l[k[i],names(which(sparse_l[i, ] > 0))])
      data=t(as.matrix(subset@assays$RNA$counts))
      if(dim(data)[1]<dim(cells_index)[1]){
        add_new=t(as.matrix(subset@assays$RNA$counts[,which(table(cells_index$A) >1)]))
        data=rbind(data,add_new)
      }
      mu_2 <- colWeightedMeans(data, w = w)
      addition=data.frame(mu_2)
      colnames(addition)=k[i]
    } else {
      addition <- data.frame(reference@assays$RNA$counts[,Cells(subset)],row.names=rownames(reference@assays$RNA$counts) )
      colnames(addition)=k[i]
    }
    merged3=cbind(merged3,addition)
  }
}

##################################


object=merged3

merged_spatial = CreateSeuratObject(counts = object, assay="Spatial")
rownames(packing_larger)=packing_larger$spot_name
coord.df = data.frame(x=packing_larger[colnames(object),"x"], y=packing_larger[colnames(object),"y"], stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) = colnames(object)

merged_spatial@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)


set.seed(19)
merged_spatial <- SCTransform(merged_spatial, assay = "Spatial", verbose = FALSE,ncells = 3000)
merged_spatial <- RunPCA(merged_spatial, assay = "SCT", verbose = FALSE) 
ElbowPlot(merged_spatial)
merged_spatial <- FindNeighbors(merged_spatial, reduction = "pca", dims = 1:50)
merged_spatial <- FindClusters(merged_spatial, verbose = FALSE,resolution = 0.5)
merged_spatial <- RunUMAP(merged_spatial, reduction = "pca", dims = 1:50)
DimPlot(merged_spatial,group.by = "seurat_clusters")

saveRDS(merged_spatial,file="/Users/naderkri/Desktop/sptype/sim/10kcells_6kspots_10um_seed15_merged2.RDS")

DimPlot(merged_spatial,group.by = "seurat_clusters")


## run ScType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ <- "/Users/naderkri/Desktop/sptype/sim/clean/reference_testing_1_20_0p4.xlsx";
gs_list <- gene_sets_prepare(db_, "Breast")
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(merged_spatial[["SCT"]])));
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(merged_spatial[["SCT"]]$scale.data) else as.matrix(merged_spatial[["SCT"]]@scale.data)
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls <- do.call("rbind", lapply(unique(merged_spatial@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(merged_spatial@meta.data[merged_spatial@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(merged_spatial@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
merged_spatial@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  merged_spatial@meta.data$sctype_classification[merged_spatial@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(merged_spatial, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')
SpatialDimPlot(merged_spatial,group.by = "sctype_classification")


DimPlot(merged_spatial, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters')








