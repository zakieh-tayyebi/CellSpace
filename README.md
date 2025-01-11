# CellSpace

CellSpace is a sequence-informed embedding method for scATAC-seq that learns a mapping of DNA *k*-mers and cells to the same space.

See our [publication](https://doi.org/10.1038/s41592-024-02274-x) for more details.

<img src="CellSpace.png" width="70%"/>

## Installation and Usage

1.  Compile the **C++** program to use as a command line tool to train a CellSpace model.

CellSpace, which uses the **C++** implementation of [StarSpace](https://github.com/facebookresearch/StarSpace) [[Wu *et al.*, 2017](https://doi.org/10.48550/arXiv.1709.03856)], builds on modern Mac OS and Linux distributions. It requires a compiler with **C++11** support and a working **make**.

Install [Boost](http://www.boost.org) library and specify the path of the library in the [makefile](cpp/makefile) (set variable **BOOST_DIR**). The default path will work if you install **Boost** by:

``` bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.63.0/source/boost_1_63_0.zip
unzip boost_1\_63_0.zip
sudo mv boost_1\_63_0 /usr/local/bin
```

Download and build CellSpace:

``` bash
git clone https://github.com/zakieh-tayyebi/CellSpace.git
cd CellSpace/cpp/
make
export PATH=$(pwd):$PATH
```

Verify that it was successfully compiled:

``` bash
CellSpace --help
```

2.  Install the **R** package to use the trained CellSpace model for downstream analysis.

Run the following commands in **R**:

``` r
install.packages("devtools")
devtools::install_github("https://github.com/zakieh-tayyebi/CellSpace.git")
library(CellSpace)
```

Installation should take only a few minutes. For details about the **R** functions, please refer to the [API](R/README.md).

3.  A tutorial on CellSpace usage can be found [here](tutorial/README.md).

## Citation

Please cite our [Nature Methods paper](https://doi.org/10.1038/s41592-024-02274-x) if you use CellSpace:

```         
Tayyebi, Z., Pine, A.R. & Leslie, C.S. Scalable and unbiased sequence-informed embedding of single-cell ATAC-seq data with CellSpace. Nature Methods 21, 1014–1022 (2024). https://doi.org/10.1038/s41592-024-02274-x
```

## Contact

-   [zakieh.tayyebi\@gmail.com](mailto:zakieh.tayyebi@gmail.com) (Zakieh Tayyebi)
-   [lesliec\@mskcc.org](mailto:lesliec@mskcc.org) (Christina S. Leslie, PhD)
