# `add_motif_db`

## Description

Computes the CellSpace embedding and activity scores of transcription factor motifs.

## Usage

``` r
add_motif_db(object, motif.db, db.name)
```

## Arguments

| Argument   | Description                                         |
|------------|-----------------------------------------------------|
| `object`   | a `CellSpace` object                                |
| `motif.db` | `PFMatrixList` or `PWMatrixList`                    |
| `db.name`  | the name of the transcription factor motif database |

## Value

a `CellSpace` object containing the motif embedding matrix, in the `motif.emb` slot, and the corresponding similarity *Z*-scores, in the `misc` slot
