# Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(dplyr)


# Import datasets
seurat_object1 <- readRDS("seurat_object1.rds")
seurat_object2 <- readRDS("seurat_object2.rds")


# Previous step if you are running the integration upon objects created from previous integrations

# Set default assay to RNA
DefaultAssay(seurat_object1) <- "RNA"
DefaultAssay(seurat_object2) <- "RNA"

# Split the dataset and create a list with one dataset per sample
seurat_object1_list <- SplitObject(seurat_object1, split.by = "orig.ident")
seurat_object2_list <- SplitObject(seurat_object2, split.by = "orig.ident")

# List of the Seurat objects
seurat_object_list <- list(seurat_object1_list, seurat_object2_list)

# SCTransform the datasets of the list
seurat_object_list <- lapply(X = seurat_object_list, FUN = SCTransform, method = "glmGamPoi",  
                             vars.to.regress = c("S.Score", "G2M.Score","Percent_mt"))


# Seurat anchor based integration

# Selection of the integration features
features <- SelectIntegrationFeatures(object.list = seurat_object_list)

# Find integration anchors
seurat_object_list <- PrepSCTIntegration(object.list = seurat_object_list,
                                         anchor.features = features) # Only when needed because SCT transformation in separate objects

anchors <- FindIntegrationAnchors(object.list = seurat_object_list, 
                                  anchor.features = features,
                                  scale = TRUE,
                                  normalization.method = "SCT",
                                  reduction = "cca") # Default, other options include "rpca", etc.

# Integration
integrated_object <- IntegrateData(anchorset = anchors)


# Necessary batch effect analysis via agrupation given the origin of samples, species... or other confounding variables

# For dataset preprocessing and clusterization change the default assay to SCT or RNA
DefaultAssay(integrated_object) <- "SCT" # or RNA

# The new dataset must be preprocessed, clusterd and annotated as the individual objects

