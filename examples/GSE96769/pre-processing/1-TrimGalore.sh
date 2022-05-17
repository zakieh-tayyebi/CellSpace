#!/bin/bash
# cutadapt
# FastQC
# TrimGalore-0.6.6

for sample in $(cat SRR_Acc_List.txt); do
  trim_galore --cores 4 --paired --fastqc --output_dir fastq/$sample $(ls fastq/$sample/*.fastq.gz)
  R1=$(ls fastq/$sample/*_R1_val_1.fq.gz); mv "$R1" "${R1%_val_1.fq.gz}_trimmed.fastq.gz"
  R2=$(ls fastq/$sample/*_R2_val_2.fq.gz); mv "$R2" "${R2%_val_2.fq.gz}_trimmed.fastq.gz"
done
