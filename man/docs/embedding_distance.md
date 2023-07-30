# `embedding_distance`

## Description

Computes distance in the embedding space based on cosine similarity.

## Usage

``` r
embedding_distance(x, y = NULL, distance = c("cosine", "angular"))
```

## Arguments

| Argument   | Description                                                                              |
|------------|------------------------------------------------------------------------------------------|
| `x`        | an embedding matrix                                                                      |
| `y`        | an embedding matrix with compatible dimensions to `x`, or `NULL`, in which case `y=x`    |
| `distance` | the distance metric, either 'cosine' or 'angular', to compute from the cosine similarity |

## Value

a matrix containing the distance between rows of `x` and `y`, computed from their cosine similarity
