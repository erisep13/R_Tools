# Monocle 2 official documentation by Trapnell Lab (https://cole-trapnell-lab.github.io/monocle-release/docs/) 

# Load packages
library(Seurat) # single-cell RNA sequencing processing
library(ggplot2) # visualization
library(sctransform)
library(glmGamPoi)
library(stringr)
library(dplyr)
library(paletteer)
library(patchwork)
library(tibble)
library(monocle) # pseudotime trajectory analysis
library(igraph) # version 1.2.11 for incompatibilities in orderCells() function


# Import Seurat processed object
seurat_object <- readRDS("path/object.rds")


# Transform Seurat object into a CellDataSet object
cds <- as.CellDataSet(seurat_object, 
                             assay = "SCT") # or RNA if not SCTransformed, expressionFamily = "negbinomial.size()" by default, recommended for most type of data


# Estimate size factors
cds <- estimateSizeFactors(cds)


# Estimate dispersions for differential expression analysis

# Patch estimateDispersionsForCellDataSet() internal function
# Copy the real function and actualize dplyr old nomenclature (group_by_() for group_by() and select_() for select())
patched_dispersions <- function(cds, modelFormulaStr, relative_expr, min_cells_detected,removeOutliers, verbose = FALSE) 
{
  if (!(("negbinomial" == cds@expressionFamily@vfamily) || 
        ("negbinomial.size" == cds@expressionFamily@vfamily))) {
    stop("Error: estimateDispersions only works, and is only needed, when you're using a CellDataSet with a negbinomial or negbinomial.size expression family")
  }
  mu <- NA
  model_terms <- unlist(lapply(str_split(modelFormulaStr, "~|\\+|\\*"), 
                               str_trim))
  model_terms <- model_terms[model_terms != ""]
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = T)
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    if (length(model_terms) > 1 || (length(model_terms) == 
                                    1 && model_terms[1] != "1")) {
      cds_pdata <- dplyr::group_by(dplyr::select(rownames_to_column(pData(cds)), 
                                                 "rowname", .dots = model_terms), .dots = model_terms)
      disp_table <- as.data.frame(cds_pdata %>% do(disp_calc_helper_NB(cds[, 
                                                                           .$rowname], cds@expressionFamily, min_cells_detected)))
    }
    else {
      cds_pdata <- dplyr::group_by(dplyr::select(rownames_to_column(pData(cds)), 
                                                 "rowname"))
      disp_table <- as.data.frame(cds_pdata %>% do(disp_calc_helper_NB(cds[, 
                                                                           .$rowname], cds@expressionFamily, min_cells_detected)))
    }
    if (!is.list(disp_table)) 
      stop("Parametric dispersion fitting failed, please set a different lowerDetectionLimit")
    disp_table <- subset(disp_table, is.na(mu) == FALSE)
    res <- parametricDispersionFit(disp_table, verbose)
    fit <- res[[1]]
    coefs <- res[[2]]
    if (removeOutliers) {
      CD <- cooks.distance(fit)
      cooksCutoff <- 4/nrow(disp_table)
      message(paste("Removing", length(CD[CD > cooksCutoff]), 
                    "outliers"))
      outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), 
                                                             names(CD)))
      res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% 
                                                  outliers == FALSE, ], verbose)
      fit <- res[[1]]
      coefs <- res[[2]]
    }
    names(coefs) <- c("asymptDisp", "extraPois")
    ans <- function(q) coefs[1] + coefs[2]/q
    attr(ans, "coefficients") <- coefs
  }
  res <- list(disp_table = disp_table, disp_func = ans)
  return(res)
}

# Internal substitution for the actual function
assignInNamespace("estimateDispersionsForCellDataSet",
                  patched_dispersions,
                  ns = "monocle")

# Load additional internal packages if they give errors
disp_calc_helper_NB <- monocle:::disp_calc_helper_NB
parametricDispersionFit <- monocle:::parametricDispersionFit

# Estimate disepersions
cds <- estimateDispersions(cds)


# Detect gene expression and filter low quality cells
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))

# Single-cell trajectory analysis

# STEP 1: DE genes for cell progress
diff_gene_test_res <- differentialGeneTest(cds[expressed_genes,],
                                           fullModelFormulaStr = "~Cluster_Labels") # variable by which DE should be made, specify "~1" for all vs all
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) # gene ids to be used in the ordering
cds <- setOrderingFilter(cds, ordering_genes) # set the ordering genes in the cds object
plot_ordering_genes(cds)

# STEP 2: Reduce data dimensionality
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree") # default reduction method

# STEP 3: Order cells along trajectory
cds <- orderCells(cds) # after doing this process for the first time we can repeat it and determine the root_state argument based on our knowledge


# Visualization
# Patching plot_cell_trajectory() function
plot_patched <- function(cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, 
                         show_backbone = TRUE, backbone_color = "black", markers = NULL, 
                         use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names = FALSE, 
                         show_state_number = FALSE, cell_size = 1.5, cell_link_size = 0.75, 
                         cell_name_size = 2, state_number_size = 2.9, show_branch_points = TRUE, 
                         theta = 0, ...) 
{
  requireNamespace("igraph")
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  lib_info_with_pseudo <- pData(cds)
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  }
  else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- reducedDimK(cds)
  }
  else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
    select(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>% 
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- minSpanningTree(cds)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_df <- dp_mst %>% igraph::as_data_frame() %>% select(source = "from", 
                                                           target = "to") %>% left_join(ica_space_df %>% select(source = "sample_name", 
                                                                                                                source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                        by = "source") %>% left_join(ica_space_df %>% select(target = "sample_name", 
                                                                                                                                             target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                                     by = "target")
  data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
    select(data_dim_1 = x, data_dim_2 = y) %>% rownames_to_column("sample_name") %>% 
    mutate(sample_state) %>% left_join(lib_info_with_pseudo %>% 
                                         rownames_to_column("sample_name"), by = "sample_name")
  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
           nrow = 2)
  }
  rot_mat <- return_rotation_mat(theta)
  cn1 <- c("data_dim_1", "data_dim_2")
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% 
                              markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData), 
      ])))
      colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
      markers_exprs <- merge(markers_exprs, markers_fData, 
                             by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                     by.y = "cell_id")
    if (use_color_gradient) {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2)) + geom_point(aes(color = value), 
                                                                      size = I(cell_size), na.rm = TRUE) + scale_color_viridis(name = paste0("value"), 
                                                                                                                               ...) + facet_wrap(~feature_label)
      }
      else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2)) + geom_point(aes(color = log10(value + 
                                                                                          0.1)), size = I(cell_size), na.rm = TRUE) + 
          scale_color_viridis(name = paste0("log10(value + 0.1)"), 
                              ...) + facet_wrap(~feature_label)
      }
    }
    else {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2, size = (value * 0.1))) + facet_wrap(~feature_label)
      }
      else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2, size = log10(value + 0.1))) + 
          facet_wrap(~feature_label)
      }
    }
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    if (use_color_gradient) {
    }
    else {
      g <- g + geom_point(aes_string(color = color_by), 
                          na.rm = TRUE)
    }
  }
  else {
    if (use_color_gradient) {
    }
    else {
      g <- g + geom_point(aes_string(color = color_by), 
                          size = I(cell_size), na.rm = TRUE)
    }
  }
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- ica_space_df %>% slice(match(mst_branch_nodes, 
                                                    sample_name)) %>% mutate(branch_point_idx = seq_len(n()))
    g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                   y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                        branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                                     size = 4, color = "white", na.rm = TRUE, branch_point_df)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_state_number) {
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }
  g <- g + monocle_theme_opts() + xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) + theme(legend.position = "top", 
                                        legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white"))
  g
}

monocle_theme_opts <- get("monocle_theme_opts", envir = asNamespace("monocle")) # Export manually from source package

# Cell-trajectory visualization
plot_patched(cds, color_by = "State") # tree-like structure plot

ggplot(as.data.frame(pData(cds)), aes(Pseudotime, reorder(Cluster_Labels, Pseudotime, median), fill=Cluster_Labels)) + geom_violin() # violin plot for cell distributions along pseudotime

seurat_object@meta.data$Pseudotime <- pData(cds)$Pseudotime # transfer pseudotimes to seurat object's metadata
FeaturePlot(seurat_object, features = "Pseudotime") + scale_color_gradientn(colors = rev(paletteer_c("grDevices::Sunset", 30))) # pseudotime over the original UMAP space

# Gene expression visualization
blast_genes <- row.names(subset(fData(cds),
                                gene_short_name %in% c("GENE1", "GENE2", "GENEn")))
plot_genes_jitter(cds[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1) # expression vs pseudotime states


expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
cds.filtered <- cds[expressed_genes,]
my_genes <- row.names(subset(fData(cds.filtered),
                             gene_short_name %in% c("GENE1", "GENE2", "GENEn")))
cds_subset <- cds.filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Cluster_Labels") # relative expression vs pseudotime coloring by eligible grouping


# Citation
citation("monocle")