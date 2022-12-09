#!/bin/bash
# samtools v1.11

cd ~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/pre-processing/var_tiles/

prepBam(){
  if [ -f aligned/$1.SAM ]; then
    samtools view -H --no-PG aligned/$1.SAM > aligned/$1\_BC.SAM
    samtools view aligned/$1.SAM | addBarcodeTag $1 >> aligned/$1\_BC.SAM
    samtools view -h -q 30 aligned/$1\_BC.SAM | \
      samtools sort -m 8G -l 9 -o aligned/$1\_BC_MQ30_posSorted.BAM -
    samtools index aligned/$1\_BC_MQ30_posSorted.BAM
    rm aligned/$1\_BC.SAM
  fi
}

g++ -std=c++11 addBarcodeTag.cpp -o addBarcodeTag
export -f prepBam
parallel -j 16 -a fastq/SRR_Acc_List.txt prepBam

ls aligned/*_BC_MQ30_posSorted.BAM > bam-files.txt
samtools merge -@ 30 -l 9 -b bam-files.txt aligned/BC_MQ30_posSorted_merged.BAM
samtools index -@ 30 aligned/BC_MQ30_posSorted_merged.BAM

