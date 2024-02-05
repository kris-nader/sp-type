
# SpType: ScType enables fast and accurate cell type identification from spatial transcriptomics data


**Article**: TBA

<p style="text-align:justify;"> <b>ScType</b> is a computational method for automated cell type annotation based on only scRNA-seq data. In this study, we adapt and showcase the application of scType, renowned for its speed, transparency, and user-friendly interface, to efficiently annotate cell types in spatial transcriptomics data.
  
Please refer to the original <a href="https://www.nature.com/articles/s41467-022-28803-w" target="_blank">ScType paper</a>  and <a href="https://github.com/IanevskiAleksandr/sc-type" target="_blank">Github</a> for more information on the method.


<br>

![alt text](https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypePlan.png)

<br>

## Quick start

```R

# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/kris-nader/sp-type/master/R/sp-type.R"); 

# Run wrapper function consists of all steps in the original sctype
# prepare gene set, calculate sctype score, add cell type annotation seurat object meta.data
sample.obj <- run_sctype(sample.obj,known_tissue_type="Brain",slot="SCT")

```

<br>

## Spatial Transcriptomics Cell type annotation example 

### Load and cluster the data


First let's load a 10x genomics dataset of sagital mouse brain slices generated using the Visium v1 chemistry. This can be done with the Seurat package. We will also go ahead and load all other recquired packages. 

```R
# Load required packages
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/kris-nader/sp-type/master/R/sp-type.R"); 

# Load demo data of sagital mouse brain anterior slice  
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
```

We will follow the recommended Seurat pipeline for processing spatial data which can be found <a href="https://satijalab.org/seurat/articles/spatial_vignette#x-visium" target="_blank">here</a>. Please note there are different protocals for analyzing spatial transcriptomics data from different platforms.

```R
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
```

### Cell Type annotation
To annotate the sample, simply apply the wrapper function to your seurat object. The function <code>run_sctype</code> uses by default the <a href="https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypeDB_full.xlsx" target="_blank">scTypeDB</a>. This consists of both positive and negative markers for cell types in the the following tissues :Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus. Users should enter the correct tissue in the parameter _known_tissue_type_. 

If your tissue of interest does not exist, feel free to use a custom marker dataset,_custom_marker_file_, which should be in the same format as the scTypeDB. Once again, we refer the reader to the original <a href="https://github.com/IanevskiAleksandr/sc-type" target="_blank">Github</a> for more information. 

<code>run_sctype</code> will output the prediction to the seurat object meta.data in _sctype_classification_.

```R
# SCT as slot/assay due to the sample being processed using the SCTransform method
brain <- run_sctype(brain,known_tissue_type="Brain",slot="SCT")
# Overlay annotation on DimPlots
b1 <- SpatialDimPlot(brain, group.by="sctype_classification")
b2 <- DimPlot(brain, group.by="sctype_classification")
b1 + b2
```


## Notes on Reproducibility
```R
sessionInfo();
[1] dplyr_1.1.4        openxlsx_4.2.5.2   HGNChelper_0.8.1   Seurat_5.0.1      
[5] SeuratObject_5.0.1 sp_2.1-1             
```

## Contact information
For any questions, please contact **Kristen Michelle Nader** (kristen.nader@helsinki.fi)

## Copyright and license

Code copyright 2024 ScType goes spatial (sp-type) , https://github.com/kris-nader/sp-type/blob/master/LICENSE
