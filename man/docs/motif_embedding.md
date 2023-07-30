# `motif_embedding`

## Description

Maps a motif to the embedding space.

## Usage

``` r
motif_embedding(object, PWM)
```

## Arguments

| Argument | Description              |
|----------|--------------------------|
| `object` | a `CellSpace` object     |
| `PWM`    | `PFMatrix` or `PWMatrix` |

## Value

a numerical vector containing the CellSpace embedding of the consensus sequence for `PWM`
