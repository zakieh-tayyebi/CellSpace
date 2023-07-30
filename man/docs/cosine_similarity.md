# `cosine_similarity`

## Description

Computes cosine similarity in the embedding space.

## Usage

``` r
cosine_similarity(x, y = NULL)
```

## Arguments

| Argument | Description                                                                           |
|----------|---------------------------------------------------------------------------------------|
| `x`      | an embedding matrix                                                                   |
| `y`      | an embedding matrix with compatible dimensions to `x`, or `NULL`, in which case `y=x` |

## Value

a matrix containing the cosine similarity between rows of `x` and `y`
