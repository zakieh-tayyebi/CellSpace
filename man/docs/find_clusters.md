# `find_clusters`

## Description

Finds clusters in a nearest neighbor graph built from the CellSpace embedding.

## Usage

``` r
find_clusters(object, graph = "cells_snn", ...)
```

## Arguments

| Argument | Description                                                                      |
|---------------------------------|---------------------------------------|
| `object` | a `CellSpace` object                                                             |
| `graph`  | name of the nearest neighbor graph in the `neighbors` slot used to find clusters |
| `...`    | arguments passed to `Seurat::FindClusters`                                       |

## Value

a `CellSpace` object with the cell clusters added to the `meta.data` slot
