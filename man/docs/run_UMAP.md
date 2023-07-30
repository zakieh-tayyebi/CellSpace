# `run_UMAP`

## Description

Computes a UMAP embedding from the CellSpace embedding.

## Usage

``` r
run_UMAP(object, emb = object@cell.emb, graph = NULL, name = "cells_UMAP", ...)
```

## Arguments

| Argument | Description                                                                                   |
|---------------------------------|---------------------------------------|
| `object` | a `CellSpace` object                                                                          |
| `emb`    | the embedding matrix used to compute the UMAP embedding                                       |
| `graph`  | name of the nearest neighbor graph in the `neighbors` slot used to compute the UMAP embedding |
| `name`   | name of the lower-dimensional embedding that will be added to the `reductions` slot           |
| `...`    | arguments passed to `Seurat::RunUMAP`                                                         |

## Value

a `CellSpace` object containing a UMAP embedding in the `reductions` slot
