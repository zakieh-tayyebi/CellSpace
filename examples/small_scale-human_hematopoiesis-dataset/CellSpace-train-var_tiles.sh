#!/bin/bash

CellSpace -thread 10 \
  -output CellSpace-results/var_tiles/CellSpace-embedding \
  -cpMat  CellSpace-inputs/GSE96769-cell_by_tile.mtx  \
  -peaks  CellSpace-inputs/GSE96769-top50K_variable_tiles.fa \
  -sampleLen 150 -ngrams 3 -exmpPerPeak 20 \
  -epoch 50 # -saveIntermediates 10

