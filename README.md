# HiCcompare

## Overview 

`HiCcompare` provides functions for joint normalization and difference detection in multiple Hi-C datasets. `HiCcompare` operates on processed Hi-C data in the form of chromosome-specific chromatin interaction matrices. `HiCcompare` is available as an R package, the major releases can be found on Bioconductor (here)[https://bioconductor.org/packages/HiCcompare/]. 

`HiCcompare` accepts three-column tab-separated text files storing chromatin interaction matrices in a sparse matrix format which are available from several sources such as the [http://aidenlab.org/data.html](http://aidenlab.org/data.html) and [http://cooler.readthedocs.io/en/latest/index.html](http://cooler.readthedocs.io/en/latest/index.html). HiCcompare is designed to give the user the ability to perform a comparative analysis on the 3-Dimensional structure of the genomes of cells in different biological states. `HiCcompare` first can jointly normalize two Hi-C datasets to remove biases between them. Then it can detect signficant differences between the datsets using a genomic distance based permutation test. The novel concept of the MD plot, based on the commonly used MA plot or Bland-Altman plot is the basis for these methods. The log **M**inus is plotted on the y axis while the genomic **D**istance is plotted on the x axis. The MD plot allows for visualization of the differences between the Hi-C datasets. 

The main functions are:
+ `hic_loess()` which performs joint `loess` normalization on the Hi-C datasets
+ `hic_diff()` which performs the difference detection process to detect significant changes between Hi-C datasets and assist in comparative analysis

Several Hi-C datasets are also included in the package.

Read the full paper describing the methods behind `HiCcompare` (here)[https://doi.org/10.1101/147850]


## Installation

First make sure you have all dependencies installed in R.

```
install.packages(c('dplyr', 'data.table', 'ggplot2', 'gridExtra', 
				   'mgcv', 'parallel', 'devtools'))

source("https://bioconductor.org/biocLite.R")
biocLite(c("InteractionSet", "GenomicRanges", "IRanges", 
		   "BiocParallel", "QDNAseq", "GenomeInfoDbData"))			   
```

To install `HiCcompare` from bioconductor open R and enter the following commands.

```
# Currently only available on the development version of Bioconductor
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("HiCcompare")
library(HiCcompare)
```


Or to install `HiCcompare` directly from the github release open R and enter the following commands.

```
library(devtools)
install_github('dozmorovlab/HiCcompare', build_vignettes = TRUE)
library(HiCcompare)
```


## Usage

First you will need to obtain some Hi-C data. Data is available from the sources listed in the overview along with many others. You will need to extract the data and read it into R as either a 3 column sparse upper triangular matrix or a 7 column BEDPE file. For more details on data extraction see the vignette included with `HiCcompare`.

Below is an example analysis using `HiCcompare`. The data in 3 column sparse upper triangular matrix format is loaded and the first step is to create a `hic.table` object using the `create.hic.table()` function. Next, the two Hi-C matrices are jointly normalized using the `hic_loess()` function. Finally, difference detection can be performed using the `hic_compare()` function. The `hic_loess()` and `hic_compare()` functions will also produce an MD plot for visualizing the differences between the datasets. 

```
# load data
library(HiCcompare)
data("HMEC.chr22")
data("NHEK.chr22")

# create the `hic.table` object
chr22.table = create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
head(chr22.table)

# Jointly normalize data for a single chromosome
hic.table = hic_loess(chr22.table, Plot = TRUE)
head(hic.table)

# input hic.table object into hic_compare
hic.table = hic_compare(hic.table, Plot = TRUE)
head(hic.table)
```

Refer to the `HiCcompare` vignette for full usage instructions. For a full explanation of the methods used in `HiCcompare` see the manuscript [here](https://doi.org/10.1101/147850).

To view the usage vignette:

`browseVignettes("HiCcompare")`

## Branches

- Master: contains the current stable release of `HiCcompare`
- supplemental: contains supplementary files and data, see Additional Vignettes section below
- test_version: contains versions of `HiCcommpare` currently in development. This version of the software may be unstable and is not reccomended for users.

## Additional Vignettes

The `HiCcompare` paper included several supplemental files that showcase some of the usage and reasoning behind the methods. Below are the titles and brief descriptions of each of these vignettes along with links to the compiled `.pdf` and the source `.Rmd` files. 

**Normalization method comparison.** 

Comparison of several Hi-C normalization techniques to display the persistence of bias in individually normalized chromatin interaction matrices, and its effect on the detection of differential chromatin interactions.

[Compiled](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S1_File.pdf)


[Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S1_File.Rmd)

**S2 File. Estimation of the IF power-law depencence.** 

Estimation of the power-law depencence between the $log_{10}-log_{10}$ interaction frequencies and distance between interacting regions. This vignette displays the reasoning behind using a power-law function for the simulation of the signal portion of Hi-C matrices.

[Compiled](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S2_File.pdf)

[Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S2_File.Rmd)

**S3 File. Estimation of the SD power-law dependence.** 

Estimation of the power-law depencence between the $log_{10}-log_{10}$ SD of interaction frequencies and distance between interacting regions. This vignette displays the reasoning behind using a power-law function for the simulation of the noise component of Hi-C matrices.

[Compiled](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S3_File.pdf)

[Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S3_File.Rmd)

**S4 File. Estimation of proportion of zeros.** 

Estimation of the depencence between the proportion of zeros and distance between interacting regions. This vignette shows distribution of zeros in real Hi-C data. The results were used for modeling the proportion of zeros in simulated Hi-C matrices with a linear function.

[Compiled](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S4_File.pdf)

[Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S4_File.Rmd)


**S5 File. Evaluation of difference detection in simulated data.** 

Extended evaluation of differential chromatin interaction detection analysis using simulated Hi-C data. Many different classifier performance measures are presented. Note: if trying to compile the source `.Rmd` this will take a long time to knit. 

[Compiled](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S5_File.pdf)

[Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S5_File.Rmd)

**S6 File. Evaluation of difference detection in real data.** 

Extended evaluation of differential chromatin interaction detection analysis using real Hi-C data. Many different classifier performance measures are presented. Note: if trying to compile the source `.Rmd` this will take a long time to knit. 

[Compiled](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S6_File.pdf)

[Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S6_File.Rmd)

**S7 File. `loess` at varying resolution.** 

Visualization of the `loess` loint normalization over varying resolutions. This vignette shows that increasing sparsity of Hi-C matrices with increasing resolution causes loess to become less useful for normalization at high resolutions. 

[Compiled](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S7_File.pdf)

[Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S7_File.Rmd)



## Citation

Please cite `HiCcompare` if you use it in your analysis.

HiCcompare: a method for joint normalization of Hi-C datasets and differential chromatin interaction detection
John Stansfield, Mikhail G. Dozmorov
bioRxiv 147850; doi: https://doi.org/10.1101/147850

## Contributions & Support

Suggestions for new features and bug reports are welcome. Please create a new [issue](https://github.com/dozmorovlab/HiCcompare/issues) for any of these or contact the author directly: [@jstansfield0](https://github.com/jstansfield0) (stansfieldjc@vcu.edu)

## Contributors

Authors: [@jstansfield0](https://github.com/jstansfield0) (stansfieldjc@vcu.edu) & [@mdozmorov](https://github.com/mdozmorov) (mikhail.dozmorov@vcuhealth.org)
