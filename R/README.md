# CellSpace (R)

This package provides tools to use and interpret a CellSpace embedding, trained with the [C++ code](../cpp/) on single-cell ATAC-seq data, in order to cluster and visualize cells, compute transcription factor activity scores, etc.

See [**Methods**](https://www.biorxiv.org/content/10.1101/2022.05.02.490310v3.full.pdf) for more details about the algorithm.

## Installing CellSpace

First, you will need to install the following **R** packages: \* [devtools](https://github.com/r-lib/devtools) \* [Seurat](https://github.com/satijalab/seurat) \* [Biostrings](https://github.com/Bioconductor/Biostrings)

Then, in order to install and use CellSpace in **R**, run the following commands:

```         
devtools::install_github("https://github.com/zakieh-tayyebi/CellSpace.git")
library(CellSpace)
```

## Documentation

Please see [the documentation](../man/) for CellSpace **R** functions. This package is not available in [CRAN](https://cran.r-project.org) yet!
