# Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(dplyr)
library(stringr)
library(biomaRt) # main tool for orthology, works as a dictionary and connects to databases (Ensembl, etc.)
library(Matrix)


# As an example: species1 = gg (Gallus gallus) & species2 = mm (Mus musculus)
# We will use Ensembl database, our mapping and alignment was made using Ensembl annotations


# Search the available datasets in Ensembl
datasets <- listDatasets(useMart("ensembl"))
datasets[grep('gallus', datasets[,1]),]
datasets[grep('mmusculus', datasets[,1]),]

# Connect to the databases
chicken <- useMart("ensembl", dataset = "ggallus_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# List chicken (species1) genes with mouse homologues
orth.gallus <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type", "mmusculus_homolog_perc_id", "mmusculus_homolog_perc_id_r1"),
  mart = chicken,
  filter = "with_mmusculus_homolog",
  values = TRUE
)

orth.gallus.1to1 <- orth.gallus[which(orth.gallus$mmusculus_homolog_orthology_type == "ortholog_one2one"),] # get only 1:1 orthologues

# List mouse (species2) genes with chicken homologues
orth.mouse <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "ggallus_homolog_ensembl_gene", "ggallus_homolog_orthology_type", "ggallus_homolog_perc_id", "ggallus_homolog_perc_id_r1"),
  mart = mouse,
  filter = "with_ggallus_homolog",
  values = TRUE
)

orth.mouse.1to1 <- orth.mouse[which(orth.mouse$ggallus_homolog_orthology_type == "ortholog_one2one"),]


# Import Seurat objects
gg <- readRDS("gg_seurat_object.rds")
mm <- readRDS("mm_seurat_object.rds")


# Filtering Seurat expression data for 1:1 orthologues in both species

# For chicken (species1)

# Extract expression matrix
gg.data <- GetAssayData(object = gg, assay = "RNA")

# Combine expression data with external gene name of Ensembl and IDs for 1 to 1 orthologues
gg.data_orth <- merge(x = orth.gallus.1to1[,1:3], y = gg.data,
                      by.x = "external_gene_name", by.y = "row.names")

# Normalize gene names and make them valid names (syntactically)
gg.data_orth$external_gene_name <- sub('[.]', '-', make.names(gg.data_orth$external_gene_name, unique=T))

# Set row names and eliminate redundant information
rownames(gg.data_orth) <- gg.data_orth$external_gene_name # set external_gene_names as rownames for the expression matrix
gg.data_orth <- gg.data_orth[,-c(1,2,3)] # eliminate the non numerical columns

# Transform the data frame into a sparse matrix
gg.data_orth <- Matrix(as.matrix(gg.data_orth), sparse = TRUE) # only numerical information can be parsed into a matrix

# Create the assay of expression for Seurat object again
gg_gene_assay <- CreateAssayObject(gg.data_orth)
gg[["RNA"]] <- gg_gene_assay

# Set the assay with the expression for 1 to 1 orthologues to work with
DefaultAssay(gg) <- "RNA"

# Repeat for mouse (species2)
mm.data <- GetAssayData(object = mm, assay = "RNA")
mm.data_orth <- merge(x = orth.mouse.1to1[,1:3], y = mm.data,
                      by.x = "external_gene_name", by.y = "row.names")
mm.data_orth$external_gene_name <- sub('[.]', '-', make.names(mm.data_orth$external_gene_name, unique=T))
rownames(mm.data_orth) <- toupper(mm.data_orth$external_gene_name) # Transform to upper case for integration with chicken - matching nomenclature
mm.data_orth <- mm.data_orth[,-c(1,2,3)]
mm.data_orth <- Matrix(as.matrix(mm.data_orth), sparse = TRUE) 
mm_gene_assay <- CreateAssayObject(mm.data_orth)
mm[["RNA"]] <- mm_gene_assay
DefaultAssay(mm) <- "RNA"

# After this, both datasets must be preprocessed and QC-filtered as done before or for the complete dataset (QC, SCTransformation, RunPCA, FindNeghbours and UMAP calculation)
# Look: seurat_scRNAseq_analysis_tool.R in https://github.com/phylobrain/R_Tools/blob/main/single-cell/seurat_scRNAseq_analysis_tool.R

# Finally, we can make a cross-species integration via anchor based integration
# Look: seurat_integration_tool.R in https://github.com/phylobrain/R_Tools/blob/main/single-cell/seurat_integration_tool.R








