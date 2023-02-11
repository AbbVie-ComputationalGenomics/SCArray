Large-scale single-cell RNA-seq data manipulation using GDS files
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

Large-scale single-cell RNA-seq data manipulation and statistical analysis with scalable implementation of generalized mixed models and principal component analysis. The package integrates the sparse matrix in Genomic Data Structure (GDS) files and the Bioconductor infrastructure framework to provide out-of-memory data storage and manipulation using the R programming language.


## Bioconductor:

v1.6.0 ([http://bioconductor.org/packages/SCArray/](http://bioconductor.org/packages/SCArray/))

Package News: [NEWS](./NEWS)


## Package Maintainer

[Xiuwen Zheng](xiuwen.zheng@abbvie.com)


## Installation

* Requires R (≥ v3.5.0), [gdsfmt](http://www.bioconductor.org/packages/gdsfmt) (≥ v1.35.4)

* Bioconductor repository
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SCArray")
```


## Examples

```R
suppressPackageStartupMessages(library(SCArray))

# the GDS file for SingleCellExperiment
fn <- system.file("extdata", "example.gds", package="SCArray")
sce <- scExperiment(fn)

sce
## class: SingleCellExperiment
## dim: 1000 850
## metadata(0):
## assays(1): counts
## rownames(1000): MRPL20 GNB1 ... RPS4Y1 CD24
## rowData names(0):
## colnames(850): 1772122_301_C02 1772122_180_E05 ... 1772122_180_B06 1772122_180_D09
## colData names(3): Cell_ID Cell_type Timepoint
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):

counts(sce)
## <1000 x 850> sparse matrix of class SC_GDSMatrix and type "double":
##              1772122_301_C02 1772122_180_E05 1772122_300_H02 ... 1772122_180_B06
##       MRPL20               3               2               3   .               0
##         GNB1              11               6              15   .               0
##        RPL22               3               5               7   .               6
##        PARK7               1               7               3   .               2
##         ENO1               8              19              20   .               7
##          ...               .               .               .   .               .
##         SSR4               0               6               3   .               5
##        RPL10              11               4               8   .               1
## SLC25A6_loc1               4               5               5   .               3
##       RPS4Y1               0               5               0   .               2
##         CD24              18               3               7   .               0
```
