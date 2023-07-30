## API

| Function                                                          | Description                                                                                     |
|-------------------------|-----------------------------------------------|
| [`CellSpace`](../man/docs/CellSpace.md)                           | Generates an object from the CellSpace class.                                                   |
| [`cosine_similarity`](../man/docs/cosine_similarity.md)           | Computes cosine similarity in the embedding space.                                              |
| [`embedding_distance`](../man/docs/embedding_distance.md)         | Computes distance in the embedding space based on cosine similarity.                            |
| [`find_neighbors`](../man/docs/find_neighbors.md)                 | Builds a nearest neighbor graph and shared nearest neighbor graph from the CellSpace embedding. |
| [`find_clusters`](../man/docs/find_clusters.md)                   | Finds clusters in a nearest neighbor graph built from the CellSpace embedding.                  |
| [`merge_small_clusters`](../man/docs/merge_small_clusters.md)     | Merges cells from small clusters with the nearest clusters.                                     |
| [`run_UMAP`](../man/docs/run_UMAP.md)                             | Computes a UMAP embedding from the CellSpace embedding.                                         |
| [`DNA_sequence_embedding`](../man/docs/DNA_sequence_embedding.md) | Maps a DNA sequence to the embedding space.                                                     |
| [`motif_embedding`](../man/docs/motif_embedding.md)               | Maps a motif to the embedding space.                                                            |
| [`add_motif_db`](../man/docs/add_motif_db.md)                     | Computes the CellSpace embedding and activity scores of transcription factor motifs.            |
