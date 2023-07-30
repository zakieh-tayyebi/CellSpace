#!/bin/bash
# TrimGalore v0.6.6

for sample in $(cat SRR_Acc_List.txt); do
  r1=$(ls fastq/$sample/*_R1.fastq.gz)
  r2=$(ls fastq/$sample/*_R2.fastq.gz)
  if [ -f $r1 ]; then
    if [ -f $r2 ]; then
      trim_galore --cores 8 --paired --fastqc --output_dir fastq/$sample $r1 $r2
      R1=$(ls fastq/$sample/*_R1_val_1.fq.gz); mv "$R1" "${R1%_val_1.fq.gz}_trimmed.fastq.gz"
      R2=$(ls fastq/$sample/*_R2_val_2.fq.gz); mv "$R2" "${R2%_val_2.fq.gz}_trimmed.fastq.gz"
    fi
  fi
done

