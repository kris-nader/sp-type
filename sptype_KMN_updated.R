library(dplyr)
library(EnhancedVolcano)
library(openxlsx)
library(HGNChelper)

sctype_source <- function(){
    # load tissue auto detect
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
    # load gene set preparation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    # load cell type annotation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    # load ScType database
    db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    return(db_)
}


run_sctype <- function(seurat_object, known_tissue_type = NULL, custom_marker_file = NULL, plot = FALSE, name = "sctype_classification",slot="RNA") {
    db_=sctype_source()
    # Check for missing arguments
    if (is.null(seurat_object)) {
        stop("Argument 'seurat_object' is missing")
    }
    if (!inherits(seurat_object, "Seurat")) {
        stop("Argument 'seurat_object' must be a Seurat object")
    }
    # Set default custom marker file
    if (is.null(custom_marker_file)) {
        custom_marker_file = db_
    }

    # Prepare gene sets
    gs_list = gene_sets_prepare(custom_marker_file, known_tissue_type)
    # Calculate scType scores
    es.max = sctype_score(scRNAseqData = seurat_object[[slot]]@scale.data,
                          scaled = TRUE,gs = gs_list$gs_positive, 
                          gs2 = gs_list$gs_negative)
    
    # Extract top cell types for each cluster
    cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
        es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
        head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
    }))
    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    seurat_object_res=seurat_object
    seurat_object_res@meta.data[name] = ""
    for(j in unique(sctype_scores$cluster)){
        cl_type = sctype_scores[sctype_scores$cluster==j,]; 
        seurat_object_res@meta.data[seurat_object_res@meta.data$seurat_clusters == j,name] = as.character(cl_type$type[1])
    }
    if(plot){
        plot_ = DimPlot(seurat_object_res, reduction = "umap", group.by = name)   
        print(plot_)
    }
    text_=paste("New metadata added: ",name)
    print(text_)
    return(seurat_object_res)
}        


get_pos_neg=function(ref_markers,cl){
    cluster_markers=ref_markers %>% filter(cluster %in% cl)
    data_EV=EnhancedVolcano(cluster_markers,rownames(cluster_markers),x ="avg_log2FC",y ="p_val_adj")
    pos_logFC=data_EV$data %>% filter(Sig=="FC_P" & avg_log2FC >0)
    neg_logFC=data_EV$data %>% filter(Sig=="FC_P" & avg_log2FC <0)
    return(list(pos_logFC,neg_logFC))
}

#' @title Read data from a zip file
#' @name sctype_spatial
#' @description saves umap and spatial umap plot pngs of annotated spatial data
#' @param data normalized seurat object of spatial 10x Visium data
#' @param ref_markers dataframe created with FindAllMarkers on normalized reference scRNA data
#' @param tissue string tissue name, capitalizede
#' @param output_name string for output folder and filenames (e.g. brca_minor)
#' @param marker_per_cluster int default: 20
#' @param neg_pos boolean, default: TRUE if filtering done on positive and negative markers separately
#' @param figure_width int default: 1000
#' @param figure_height int default: 550
#' @param output_folder name of the output file for plots default: "figures"
#' @return 
#' @export
#' @examples
retrieve_markers=function(data, ref_markers, tissue, marker_per_cluster = NULL){
    
    ScTypeDB = data.frame(tissueType = rep(tissue, length(unique(ref_markers$cluster))),
                           cellName = unique(ref_markers$cluster),
                           geneSymbolmore1 = "",
                           geneSymbolmore2 = "")
    
    rownames(ScTypeDB)=ScTypeDB$cellName
    for (cl in as.character(unique(ref_markers$cluster))) {
        df_=get_pos_neg(ref_markers,cl)
        df_pos=df_[[1]][order(df_[[1]]$p_val_adj,(df_[[1]]$avg_log2FC) * -1),] 
        df_neg=df_[[2]][order(df_[[2]]$p_val_adj,(df_[[2]]$avg_log2FC)),] 
        if (is.null(marker_per_cluster)) {
          positive_markers=paste0(df_pos$gene,collapse = ",")
          negative_markers=paste0(df_neg$gene,collapse = ",")
        }
        else {
          positive_markers=paste0(df_pos$gene[1:marker_per_cluster],collapse = ",")
          negative_markers=paste0(df_neg$gene[1:marker_per_cluster],collapse = ",")
        }
        ScTypeDB[cl,"geneSymbolmore1"]=positive_markers
        ScTypeDB[cl,"geneSymbolmore2"]=negative_markers
        }
    rownames(ScTypeDB)=NULL
    pwd_=getwd()
    dir.create(paste0(pwd_,"/temp"))
    excel_filename <- paste0(pwd_,"/temp/","sctypeDB_", marker_per_cluster, ".xlsx")
    write.xlsx(ScTypeDB, file = excel_filename)
}

run_plot_sctype=function(db_=NULL,tissue,data,figure_width=1000,figure_height=550,
                         output_name = paste0(tissue, "_sptype_output"),
                         output_folder = "figures", marker_per_cluster = NULL, saveRDS = TRUE,
                         pt.size.factor = 2.4, st_method = "visium"){
    data_annotated=run_sctype(data,known_tissue_type=tissue,custom_marker_file=db_,name="customclassif",slot="SCT")
    Idents(data_annotated) <- data_annotated@meta.data$customclassif
    if (st_method == "visium") {
      p1 <- DimPlot(data_annotated, reduction = "umap", label = FALSE, repel = TRUE, group.by = 'customclassif') 
      p2 <- SpatialDimPlot(data_annotated, group.by = 'customclassif', pt.size.factor = pt.size.factor, alpha = c(1, 1), label = FALSE)
    }
    else {
      p1 <- DimPlot(data_annotated, reduction = "umap", label = FALSE, repel = TRUE, group.by = 'customclassif', raster = FALSE) 
      p2 <- ImageDimPlot(data_annotated, group.by = 'customclassif')
    }
    #save plots
    dir.create(output_folder, recursive = TRUE)
    png(file = paste0("./", output_folder, "/", output_name, "_sptype_",
                      marker_per_cluster, ".png"),width=figure_width, height=figure_height)
    print(p1 + p2)
    dev.off()
    if (saveRDS) {
      if (!is.null(marker_per_cluster)) {
        rds_filename <- paste0(output_name, "_", marker_per_cluster, ".RDS")
      }
      else {
        rds_filename <- paste0(output_name, ".RDS")
      }
      saveRDS(data_annotated,file = rds_filename)
    }
    return(data_annotated)
}


#load("brca_her2_ref_markers_minor_normalized.rds") #or major
#load("brca.rds")

#retrieve_markers(data, ref_markers, tissue="Breast")

#load("/Users/naderkri/Downloads/brca_1.rds")
#SpatialDimPlot(data, group.by = 'patho_annot', pt.size.factor = 2.4, alpha = c(1, 1), label = FALSE)
