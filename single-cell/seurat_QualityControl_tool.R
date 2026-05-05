# Tool to visualize and establish the "optimal" thresholds for QC
# QC thresholds depend on the species, the tissue and the sequencing method used to obtain the data
# Single-nuclei data tend to have less count of nFeatures and nCounts than single-cell data, we also expect to see less %MT

# Load packages
library(Seurat)
library(ggplot2)

# Import data / Construct the Seurat object
seurat_object <- readRDS("seurat_object.rds")

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
# The first and most important step for QC threshold election is watching the cell distribution

Idents(seurat_object) <- "Sample1"
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "Percent_mt"), group.by = "orig.ident")


FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "Percent_mt", group.by = "orig.ident")

# We expect to see a linear correlation between the number of genes (nFeature_RNA) and the number of UMIs (nCounts_RNA)
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")

# Other types of distribution visualizations
seurat_object@meta.data %>% 
  ggplot(aes(color=Sample_Name, x=nFeature_RNA, fill= Sample_Name)) + # or nCounts_RNA / Percent_mt
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")

seurat_object@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=Percent_mt)) + # Linear relationship between features and UMIs and %MT expression
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic()

# Quantile election for our data and quantile filtering (2%-98%, for example)
nFeature_low <- quantile(seurat_object$nFeature_RNA, 0.02) # also for nCount_RNA or Percent_mt
nFeature_high <- quantile(seurat_object$nFeature_RNA, 0.98)

# QC filtering
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 600 & nFeature_RNA < 3500 & Percent_mt < 0.5) # We can decide manually the thresholds or include the percentile calculated thresholds


# Next we continue with the preprocessing of the data