# R_Tools

This repository contains tools for working with R (v4.5.2) in RStudio (v0.18.0), running under Ubuntu 24.04.4 LTS. 

## Contents

- **single-cell**: Contains tools for single-cell RNA sequencing analysis, mainly based on [Seurat](https://satijalab.org/seurat/) pipeline.
  - *monocle2_pseudotime_tool.R*: Pseudotime trajectory analysis for single-cell RNA sequencing data previously processed and annotated with Seurat using Monocle 2 package version (v2.38.0).
  - *seurat_biomart_orthology_tool.R*: One to one orthologue extraction and filtering for cross-species comparations using the `biomaRt`package (v2.66.2) with [Ensembl](https://www.ensembl.org/index.html) gene annotations.
  - *seurat_integration_tool.R*: Dataset integration using Seurat anchor based integration method.
  - *seurat_scRNAseq_analysis_tool.R*: Basic pipeline for Seurat object construction from a count matrix (obtained via Fluent Illumina technology), quality control, preprocessing and clustering, with basic visualization. Preprocessing, scaling and normalization using the SCTransform pipeline (v0.4.3) and Seurat (v5.4.0) base functions.
  
  

## Session and package information

### Specific analysis packages

- `Seurat` (v5.4.0): [Hao, Y., *et al.*, (2024)](https://www.nature.com/articles/s41587-023-01767-y#Abs1)
- `monocle` (v2.38.0): [Trapnell, C., *et al.* (2014)](https://www.nature.com/articles/nbt.2859) & [Qiu, X., *et al.*, (2017)](https://www.nature.com/articles/nmeth.4402)
- `sctranform` (v0.4.3): [Hafemeister, C. & Satija, R., (2019)](https://doi.org/10.1186/s13059-019-1874-1), [Choudhary, S. & Satija, R., (2022)](https://doi.org/10.1186/s13059-021-02584-9) & [CRAN package](https://doi.org/10.32614/CRAN.package.sctransform)
- `glmGamPoi` (v1.22.0): [Ahlmann-Eltze, C. & Huber, W., (2020)](https://doi.org/10.1093/bioinformatics/btaa1009) & [Bioconductor Package](https://doi.org/doi:10.18129/B9.bioc.glmGamPoi)
- `biomaRt` (v2.66.2): [Durnick, S., *et al.*, (2005)](https://doi.org/10.1093/bioinformatics/bti525), [Durnick, S., *et al.*, (2009)](https://doi.org/10.1038/nprot.2009.97) & [Bioconductor Package](https://doi.org/doi:10.18129/B9.bioc.biomaRt)

### General R language packages

- `ggplot2` (v4.0.2): [Wickham, H., (2016)](https://doi.org/10.1007/978-3-319-24277-4), [ggplot2 tidyverse](https://ggplot2.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.ggplot2)
- `dplyr` (v1.2.0): [Wickham, H., *et al.*, (2026)](https://dplyr.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.dplyr)
- `stringr` (v1.6.0): [Wickham, H., (2025)](https://stringr.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.stringr)
- `tibble` (v1.6.0): [Müller, K. & Wickham, H. (2025)](https://tibble.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.tibble)
- `patchwork`(v1.3.2): [Pedersen, T. (2025)](https://patchwork.data-imaginist.com/) & [CRAN package](https://doi.org/10.32614/CRAN.package.patchwork)
- `paletteer` (v1.7.0): [Hvitfeldt, E. (2021)](https://github.com/EmilHvitfeldt/paletteer), [paletteer official website](https://emilhvitfeldt.github.io/paletteer/index.html) & [CRAN package](https://doi.org/10.32614/CRAN.package.paletteer)
