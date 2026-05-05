# Single-cell RNA sequencing analysis in Seurat starting from count matrix

# Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(dplyr)
library(stringr)


# Data import and construction of the count matrix
counts <- ReadMtx(mtx = "matrix.mtx", cells = "barcodes.tsv", features = "features.tsv") # creation of a sparse matrix with cells, features and counts

# Create Seurat objects
seurat_object <- CreateSeuratObject(counts=counts) # more arguments available

# Addition of metadata
seurat_object <- AddMetaData(seurat_object, "Sample1", col.name = "Sample_Name") # add as many meta data as necessary to describe your data


# Quality Control

# Mitochondrial percentage calculation
seurat_object[["Percent_mt"]] <- PercentageFeatureSet(seurat_object, pattern = "mt-") # for mouse
seurat_object[["Percent_mt"]] <- PercentageFeatureSet(seurat_object, pattern = "MT-") # for human
seurat_object[["Percent_mt"]] <- PercentageFeatureSet(seurat_object, features = c("COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB", "ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6", 
                                                                   "ENSGALG00010000002", "ENSGALG00010000003", "ENSGALG00010000004", "ENSGALG00010000005", "ENSGALG00010000006", 
                                                                   "ENSGALG00010000008", "ENSGALG00010000009", "ENSGALG000100000010", "ENSGALG000100000012", "ENSGALG000100000013", 
                                                                   "ENSGALG000100000014", "ENSGALG000100000015", "ENSGALG000100000016", "ENSGALG000100000018", "ENSGALG000100000019", 
                                                                   "ENSGALG000100000021", "ENSGALG000100000025", "ENSGALG000100000027", "ENSGALG000100000030", "ENSGALG000100000031", 
                                                                   "ENSGALG000100000032", "ENSGALG000100000035", "ENSGALG000100000036", "ENSGALG000100000038")) # for chicken

# Distribution visualization
Idents(seurat_object) <- "Sample1"
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "Percent_mt"), group.by = "orig.ident")
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "Percent_mt", group.by = "orig.ident")
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") # closer to 1 better

# QC filtering
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & Percent_mt < 5) # thresholds are changed based on the samples, this are some standards

# SCTransoform without out-regression
seurat_object <- SCTransform(seurat_object, verbose = FALSE)

# Cell cycle scoring
fgm.s.genes <- c("GINS1", "RFC1", "WXO1", "ENSGALG00000011747", "H2AFJ", "ENSGALG00000042491", "HIST1H110") # FGM and RSG's phase S gene list
fgm.g2m.genes <- c("BUB1B", "NCAPG", "NCAPH", "CDC25B", "SMC2", "SPC25") # FGM and RSG's phase G2m gene list
s.genes <- toupper(c(cc.genes$s.genes)) # Convert to upper or title-like case depending on the gene nomenclature of the species
g2m.genes <- toupper(c(cc.genes$g2m.genes))

seurat_object <- CellCycleScoring(seurat_object, s.features = c(s.genes, fgm.s.genes), g2m.features = c(g2m.genes, fgm.g2m.genes), set.ident = TRUE, nbin = 10)

# SCTransform with out-regression if needed/wanted
seurat_object <- SCTransform(seurat_object, vars.to.regress = c("S.Score", "G2M.Score", "Percent_mt"), verbose = FALSE)

# Preprocessing and clustering
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object, ndims=50) # Choose optimal PCs/dimensions for UMAP calculation
seurat_object <- FindNeighbors(seurat_object, verbose = FALSE, dims = 1:30) # Set the dimensions chosen in PCA's Elbow plot
seurat_object <- FindClusters(seurat_object, resolution = c(0.02,0.5,1,1.5,2)) # Try different resolutions for better model fitting

set.seed(123) # Set random seed for reproducibility of UMAP plot due to randomness
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:30) # Set the dimensions chosen in PCA's Elbow plot, consistent with number of PCs in FindNeighbors()


# Visualization

# Visualization of cell cycle scoring and phase inferring
DimPlot(seurat_object, group.by = "Phase") # Inferred cell cycle phase during CellCycleScoring() function
FeaturePlot(seurat_object, features = c("S.Score", "G2M.Score", "Percent_mt"))
DotPlot(seurat_object, features =  c("S.Score", "G2M.Score", "Percent_mt"))
VlnPlot(seurat_object, features =  c("S.Score", "G2M.Score", "Percent_mt"))

# Clustering visualization
DimPlot(seurat_object, group.by = "SCT_snn_res.0.02", label = T, label.size = 4, repel = T) 


# Annotations

# Find variable markers among clusters
# Preprocessing of SCT assay
seurat_object <- PrepSCTFindMarkers(seurat_object, assay = "SCT", verbose = TRUE) # Optional, when normalizing with SCT and after integration itt may be necessary

markers <- FindAllMarkers(object = seurat_object,
                          only.pos = TRUE,
                          verbose = FALSE,
                          group.by = "SCT_snn_res.1", # Set the resolution you are going to work with
                          min.pct = 0.25,
                          logfc.threshold = 0.5) # Threshold for better marker finding
# Top n genes for each cluster
topn_genes <- markers %>% group_by(cluster) %>% top_n(n, avg_log2FC)
topn_genes <- topn_genes %>%
  filter(!grepl("^ENSGAL", gene)) # For filtering of genes without a gene name (optional)

# Top genes visualization
DoHeatmap(seurat_object, group.by = "SCT_snn_res.1", features=topn_genes$gene, raster=F, label=F) +
  theme(axis.text.y = element_text(size = 8))
DotPlot(seurat_object, features = unique(topn_genes$gene), group.by = "SCT_snn_res.1") + coord_flip() +
  theme(axis.text.y = element_text(size = 8))

# Rename Idents
Idents(seurat_object) <- "SCT_snn_res.1" # Resolution you are working with
seurat_object <- RenameIdents(seurat_object,
                                 c("0" = "Name1",
                                   "1" = "Name2",
                                   "2" = "Name3",
                                   "3" = "Name2",
                                   "4" = "Name1",
                                   "5" = "Name2",
                                   "6" = "Name4")) # Match cluster number with the new cluster name
