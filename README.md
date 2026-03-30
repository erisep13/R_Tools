# R_Tools

This repository contains tools for working with R (v4.5.2) in RStudio (v0.18.0), running under Ubuntu 24.04.4 LTS. 

## Contents

- **single-cell**: Contains tools for single-cell RNA sequencing analysis, mainly based on [Seurat](https://satijalab.org/seurat/) pipeline.
  - *monocle2_pseudotime_tool.R*: Pseudotime trajectory analysis for single-cell RNA sequencing data previously processed and annotated with Seurat using Monocle 2 package version (v2.38.0).

## Session and package information

### Specific analysis packages

- `Seurat` (v5.4.0): [Hao, Y., *et al.*, (2024)](https://www.nature.com/articles/s41587-023-01767-y#Abs1)
- `monocle` (v2.38.0): [Trapnell, C., *et al.* (2014)](https://www.nature.com/articles/nbt.2859) & [Qiu, X., *et al.*, (2017)](https://www.nature.com/articles/nmeth.4402)
- `sctranform` (v0.4.3): [Hafemeister, C. & Satija, R., (2019)](https://doi.org/10.1186/s13059-019-1874-1), [Choudhary, S. & Satija, R., (2022)](https://doi.org/10.1186/s13059-021-02584-9) & [CRAN package](https://doi.org/10.32614/CRAN.package.sctransform)
- `glmGamPoi` (v1.22.0): [Ahlmann-Eltze, C. & Huber, W., (2020)](https://doi.org/10.1093/bioinformatics/btaa1009) & [Bioconductor Package](https://doi.org/doi:10.18129/B9.bioc.glmGamPoi) 

### General R language packages

- `ggplot2` (v4.0.2): [Package tidyverse website](https://ggplot2.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.ggplot2)
- `dplyr` (v1.2.0): [Package tidyverse website](https://dplyr.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.dplyr)
- `stringr` (v1.6.0): [Package tidyverse website](https://stringr.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.stringr)
- `tibble` (v1.6.0): [Package tidyverse website](https://tibble.tidyverse.org/) & [CRAN package](https://doi.org/10.32614/CRAN.package.tibble)
- `patchwork`(v1.3.2): [Official package website](https://patchwork.data-imaginist.com/) & [CRAN package](https://doi.org/10.32614/CRAN.package.patchwork)
- `paletteer` (v1.7.0): [CRAN package](https://doi.org/10.32614/CRAN.package.paletteer)
