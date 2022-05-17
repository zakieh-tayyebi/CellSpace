#!/bin/bash
# samtools-1.11

prepBam(){
  aligned=aligned/$sample.SAM
  if [ -f $aligned ]; then
    samtools view -H --no-PG aligned/$sample.SAM > aligned/$sample\_BC.SAM
    samtools view aligned/$sample.SAM | ./addBarcodeTag $sample >> aligned/$sample\_BC.SAM
    samtools view -h -q 30 aligned/$sample\_BC.SAM | \
      samtools sort -m 8G -l 9 -o aligned/$sample\_BC_MQ30_posSorted.BAM -
    samtools index aligned/$sample\_BC_MQ30_posSorted.BAM
    rm aligned/$sample\_BC.SAM
  fi
}
export -f prepBam
parallel -a SRR_Acc_List.txt prepBam

ls aligned/*_BC_MQ30_posSorted.BAM > bam-files.txt
samtools merge -@ 30 -l 9 -b bam-files.txt aligned/BC_MQ30_posSorted_merged.BAM
samtools index -@ 30 aligned/BC_MQ30_posSorted_merged.BAM
rm bam-files.txt aligned/*_BC_MQ30_posSorted.BAM*
