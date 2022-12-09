#!/bin/bash

CellSpace -thread 10 \
  -output CellSpace-results/all_peaks/CellSpace-embedding \
  -cpMat  CellSpace-inputs/GSE96769-cell_by_peak.mtx  \
  -peaks  CellSpace-inputs/GSE96769-filtered_peaks.fa \
  -sampleLen 150 -ngrams 3 -exmpPerPeak 20 \
  -epoch 50 # -saveIntermediates 10

