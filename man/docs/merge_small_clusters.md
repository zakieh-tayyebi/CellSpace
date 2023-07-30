# `merge_small_clusters`

## Description

Merges cells from small clusters with the nearest clusters.

## Usage

``` r
merge_small_clusters(
  object,
  clusters,
  min.cells = 10,
  graph = "cells_snn",
  seed = 1
)
```

## Arguments

| Argument    | Description                                                                                                      |
|---------------------------------|---------------------------------------|
| `object`    | a `CellSpace` object                                                                                             |
| `clusters`  | a vector of cluster labels, or the name of a column in the `meta.data` slot containing cluster labels            |
| `min.cells` | any cluster with fewer cells than `min.cells` will be merged with the nearest cluster                            |
| `graph`     | a nearest neighbor graph, or the name of a nearest neighbor graph in the `neighbors` slot, used to find clusters |

## Value

new cluster labels
