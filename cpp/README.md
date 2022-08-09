# CellSpace (C++)

This code is an extension of [StarSpace](https://github.com/facebookresearch/StarSpace) (see [MIT license](StarSpace-MIT-License.md)) and incorporates the CellSpace algorithm for scATAC-seq data.

See [the methods](https://www.biorxiv.org/content/early/2022/05/20/2022.05.02.490310.full.pdf) for more details about the algorithm.

## Building CellSpace

CellSpace builds on modern Mac OS and Linux distributions. It requires a compiler with **C++11** support and a working **make**.

You need to install <a href=http://www.boost.org/>Boost</a> library and specify the path of the library in the makefile (**BOOST_DIR**).

In order to build and use CellSpace, run the following commands:

    git clone https://github.com/zakieh-tayyebi/CellSpace.git
    cd CellSpace/cpp/
    make
    export PATH=$(pwd):$PATH

## Input Parameters

    CellSpace ...
    
    The following arguments are mandatory:
      -output             prefix of the output
      -cpMat              sparse cell by peak/tile count matrix (.mtx)
      -peaks              multi-fasta file containing peak/tile DNA sequences with the order they appear in the corresponding count matrix (.fa)
    
    The following arguments are optional:
      -dim                size of embedding vectors [default=30]
      -ngrams             max length of k-mer ngram [default=3]
      -k                  k-mer length [default=8]
      -sampleLen          length of the sequences randomly sampled from the peak/tile DNA sequences (integer or 'given') [default=150]
      -exmpPerPeak        number of training examples per peak/tile [default=20]
      -epoch              number of epochs [default=20]
      -margin             margin parameter in hinge loss [default=0.05]
      -bucket             number of buckets [default=2000000]
      -label              labels prefix [default='__label__']
      -lr                 learning rate [default=0.01]
      -maxTrainTime       max train time (seconds) [default=8640000]
      -negSearchLimit     number of negative labels sampled per dataset [default=50]
      -maxNegSamples      max number of negatives in a batch update [default=10]
      -p                  the embedding of an entity equals the sum of its M feature embedding vectors devided by M^p [default=0.5]
      -initRandSd         initial values of embeddings are randomly generated from normal distribution with mean=0 and standard deviation=initRandSd [default=0.001]
      -batchSize          size of mini batch in training [default=5]
      -saveIntermediates  save intermediate models or only the final model (integer or 'final') [default='final']
      -thread             number of threads [default=10]

* In the count matrix provided by `-cpMat`, rows are associated with cells and columns with accessibility events (peaks or tiles). Events that are not accessible in any cells (i.e. all-zero columns) must be excluded. In **R**, a sparse matrix of class **dgCMatrix** can be written into a **.mtx** file with **Matrix::writeMM**.
* In order to train CellSpace on multiple datasets, their input files must be provided as a list: `-cpMat data1-counts.mtx data2-counts.mtx -peaks data1-peaks.fa data2-peaks.fa`
* If `-sampleLen given` is specified, the entire peak/tile sequence is used.
