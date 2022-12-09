#!/bin/bash
# sratoolkit v2.11.0-ubuntu64

cd ~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/pre-processing/var_tiles/

i=0
for sample in $(cat fastq/SRR_Acc_List.txt); do
  let i++; echo ==================== S$i $sample >> fastq/fastq-dump.std
  fastq-dump --split-3 --dumpbase --skip-technical --clip --read-filter pass \
             --outdir fastq $sample >> fastq/fastq-dump.std 2>&1
  gzip fastq/*.fastq

  R1=fastq/$sample\_pass_1.fastq.gz
  R2=fastq/$sample\_pass_2.fastq.gz
  if [ -f $R1 ]; then
    if [ -f $R2 ]; then
      mkdir fastq/$sample
      mv $R1 fastq/$sample/$sample\_S$i\_R1.fastq.gz
      mv $R2 fastq/$sample/$sample\_S$i\_R2.fastq.gz
    fi
  fi
done
rm fastq/*_pass_*.fastq.gz

