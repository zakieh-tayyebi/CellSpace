# `find_neighbors`

## Description

Builds a nearest neighbor graph and shared nearest neighbor graph from the CellSpace embedding.

## Usage

``` r
find_neighbors(
  object,
  n.neighbors = 30,
  emb = object@cell.emb,
  emb.name = "cells",
  ...
)
```

## Arguments

| Argument      | Description                                                           |
|---------------|-----------------------------------------------------------------------|
| `object`      | a `CellSpace` object                                                  |
| `n.neighbors` | the number of nearest neighbors for the KNN algorithm                 |
| `emb`         | the embedding matrix used to create the nearest neighbor graphs       |
| `emb.name`    | prefix for the graph names that will be added to the `neighbors` slot |
| `...`         | arguments passed to `Seurat::FindNeighbors`                           |

## Value

a `CellSpace` object containing nearest neighbor and shared nearest neighbor graphs in the `neighbors` slot
