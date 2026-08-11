# Load packages
library(clustree) # https://doi.org/10.1093/gigascience/giy083

# Load a Seurat object
seurat_object <- readRDS("seurat_object.rds")

# If resolutions previously calculated, better clean those columns
cols <- grep("^SCT_snn_res", colnames(seurat_object@meta.data), value = TRUE)
seurat_object@meta.data[, cols] <- NULL

# Cluster again selecting increasing numbers of resolutions you want to analyse
# Calculate several resolutions
seurat_object <- RunPCA(seurat_object)
seurat_object <- FindNeighbors(seurat_object, verbose = FALSE, dims = 1:n)
seurat_object <- FindClusters(seurat_object, 
                              resolution = c(0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25)) # select the resolutions you are checking
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:n)

# Build a clustering tree
clustree(seurat_object, 
         prefix = "SCT_snn_res.") # Prefix of the resolution, everything before the number of res


# See the stability of the clusters
clustree(seurat_object, prefix = "SCT_snn_res.", node_colour = "sc3_stability")

# Color by other variables in the metadata, for example Percent_mt
clustree(seurat_object, prefix = "SCT_snn_res.", node_colour = "Percent_mt", node_colour_aggr = "mean")

# Color by gene expression
clustree(seurat_object, prefix = "SCT_snn_res.",
         node_colour = "CORIN", node_colour_aggr = "median") # here we color by CORIN gene median value

# Overlying clustering trees
clustree_overlay(seurat_object, prefix = "SCT_snn_res.", red_dim = "pca", 
                 x_value = "pca_1", y_value = "pca_2") # also over umap dimensions

# Citing clustree
citation("clustree")
