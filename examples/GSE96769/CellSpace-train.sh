#!/bin/bash

CellSpace \
  -output CellSpace-results/GSE96769-embedding.tsv \
  -cpMat CellSpace-inputs/cell-by-tile-bin.mtx  \
  -peaks CellSpace-inputs/top50K-variable-tiles.fa \
  -sampleLen 150 -ngrams 3 -exmpPerPeak 20 -epoch 20

