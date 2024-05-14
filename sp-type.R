#########################################################
##         Get functions for sctype goes spatial       ##
#########################################################
#
# GNU General Public License v3.0 (https://github.com/kris-nader/sp-type/blob/main/LICENSE)
#
# Written by Kristen Michelle Nader <kristen.nader@helsinki.fi> February 2024
#
# Functions on this page:
# sctype_source,run_scType




#' @title sctype source files
#' @name sctype_source
#' @description loads sctype functions needed for an automated cell type annotation . 
#' @details none
#' @param none 
#' @return original ScType database
#' @export
#' @examples
#' db_=sctype_source()
#' 


sctype_source <- function(){
    # load gene set preparation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    # load cell type annotation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    # load ScType database
    db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    return(db_)
}

#' @title Run sctype analysis on Seurat object
#' @name run_scType
#' @description run an automated cell type annotation on spatial data-- spot type annotation
#' @details Useful to get an idea of different cells in the sample
#' @param seurat_object A Spatial Seurat object
#' @param known_tissue_type The tissue type of the input data should match what is in referenceDB
#' @param custom_marker_file Path to the custom marker file (optional)
#' @param plot Whether to plot the results (default is FALSE)
#' @param name The name of the metadata column to store the scType results (default is "sctype_classification")
#' @return A modified copy of the input Seurat object with a new metadata column
#' 
#' @import sctype source code
#' @import Seurat DimPlot
#' 
#' @examples
#' seurat_object=run_scType(seurat_object,known_tissue_type = "Hippo",slot = "SCT",name = "sctype-goes-spatial",custom_marker_file = "./ref_markers_brain.xlsx")
#' 
#' @export
#'

run_sctype <- function(seurat_object, known_tissue_type = NULL, custom_marker_file = NULL, plot = FALSE, name = "sctype_classification",slot="SCT") {
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
    es.max = sctype_score(scRNAseqData = seurat_object[[slot]]$scale.data,
                          scaled = TRUE,gs = gs_list$gs_positive, 
                          gs2 = gs_list$gs_negative)
    
    # Extract top cell types for each cluster
    cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
        es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
        head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
    }))
    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/40] = "Unknown"
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


get_pos_neg=function(ref_markers, cl, min_pct_diff){
    cluster_markers=ref_markers %>% filter(cluster %in% cl)
    cluster_markers <- cluster_markers[abs(cluster_markers$pct.1 - cluster_markers$pct.2) > min_pct_diff,]
    data_EV=EnhancedVolcano(cluster_markers,rownames(cluster_markers),x ="avg_log2FC",y ="p_val_adj")
    pos_logFC=data_EV$data %>% filter(Sig=="FC_P" & avg_log2FC >0)
    neg_logFC=data_EV$data %>% filter(Sig=="FC_P" & avg_log2FC <0)
    return(list(pos_logFC,neg_logFC))
}


retrieve_markers=function(ref_markers, tissue, marker_per_cluster = NULL, min_pct_diff = 0.3){
    
    ScTypeDB = data.frame(tissueType = rep(tissue, length(unique(ref_markers$cluster))),
                           cellName = unique(ref_markers$cluster),
                           geneSymbolmore1 = "",
                           geneSymbolmore2 = "")
    
    rownames(ScTypeDB)=ScTypeDB$cellName
    for (cl in as.character(unique(ref_markers$cluster))) {
        df_=get_pos_neg(ref_markers,cl, min_pct_diff = min_pct_diff)
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
    
    return(ScTypeDB)
}
