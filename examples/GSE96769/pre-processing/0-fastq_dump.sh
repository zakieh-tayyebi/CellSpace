#!/bin/bash
# sratoolkit.2.11.0-ubuntu64

mkdir fastq
i=0
for sample in $(cat SRR_Acc_List.txt); do # SRA accession list from GEO
  let i++

  fastq-dump --split-3 --dumpbase --skip-technical --clip --read-filter pass --outdir fastq $sample
  gzip fastq/*.fastq

  R1=fastq/$sample\_pass_1.fastq.gz; R2=fastq/$sample\_pass_2.fastq.gz
  if [ -f $R1 ]; then if [ -f $R2 ]; then
    mkdir fastq/$sample
    mv $R1 fastq/$sample/$sample\_S$i\_R1.fastq.gz
    mv $R2 fastq/$sample/$sample\_S$i\_R2.fastq.gz
  fi; fi
done
