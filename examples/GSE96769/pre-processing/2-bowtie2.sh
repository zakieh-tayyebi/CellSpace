#!/bin/bash
# bowtie2-2.4.2-linux-x86_64

mkdir aligned
for sample in $(cat SRR_Acc_List.txt); do
  R1=$(ls fastq/$sample/*_R1_trimmed.fastq.gz)
  R2=$(ls fastq/$sample/*_R2_trimmed.fastq.gz)
  if [ -f $R1 ]; then if [ -f $R2 ]; then
    bowtie2 --threads 8 --mm --maxins 2000 --very-sensitive --no-unal --no-mixed \
      -x Bowtie2_index/hg19 -1 $R1 -2 $R2 -S aligned/$sample.SAM
  fi; fi
done
