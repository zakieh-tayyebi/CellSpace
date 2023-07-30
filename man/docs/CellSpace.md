# `CellSpace`

## Description

Generates an object from the `CellSpace` class.

## Usage

``` r
CellSpace(
  emb.file,
  cell.names = NULL,
  meta.data = NULL,
  project = NULL,
  similarity = "cosine",
  p = 0.5,
  label = "__label__"
)
```

## Arguments

| Argument     | Description                                                                                     |
|-----------------|-------------------------------------------------------|
| `emb.file`   | the .tsv output of CellSpace containing the embedding matrix for cells and *k*-mers             |
| `cell.names` | vector of unique cell names                                                                     |
| `meta.data`  | a `data.frame` containing meta-information about each cell                                      |
| `project`    | title of the project                                                                            |
| `similarity` | the similarity function in hinge loss                                                           |
| `p`          | the embedding of an entity equals the sum of its `M` feature embedding vectors divided by `M^p` |
| `label`      | cell label prefix                                                                               |

## Value

a new `CellSpace` object
