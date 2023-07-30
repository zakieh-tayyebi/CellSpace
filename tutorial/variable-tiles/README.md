Pre-processing the 'Small-scale human hematopoiesis dataset':

1.  [Downloading FASTQ files](1-fastq-dump.sh),
2.  [Trimming adapter sequences](2-TrimGalore.sh),
3.  [Aligning reads to the reference genome](3-bowtie2.sh),
4.  [Filtering and sorting the alignments](4-samtools.sh),
5.  [Filtering fragments and cells with ArchR (quality-control and doublet removal)](5-ArchR.R).

The result is an ArchR object containing the pre-processed scATAC-seq data ('archr.obj').
