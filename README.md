# [HiCcompare](https://dozmorovlab.github.io/HiCcompare/)

<!-- badges: start -->
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/HiCcompare.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/HiCcompare)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
<!-- badges: end -->

## Overview 

`HiCcompare` provides functions for joint normalization and difference detection in multiple Hi-C datasets. `HiCcompare` operates on processed Hi-C data in the form of chromosome-specific chromatin interaction matrices. [HiCcompare is available on Bioconductor](https://bioconductor.org/packages/HiCcompare/). 

If you have more than two Hi-C datasets, please see our other package, [multiHiCcompare](https://bioconductor.org/packages/multiHiCcompare/).

`HiCcompare` accepts three-column tab-separated text files storing chromatin interaction matrices in a sparse matrix format which are available from several sources such as the [Aiden Lab](http://aidenlab.org/data.html) (`.hic` files) and the [Mirnylab FTP site](http://cooler.readthedocs.io/en/latest/index.html) (`.cool` files). `HiCcompare` performs differential chromatin interaction analysis between two biological conditions, one Hi-C matrix per condition. 

First, `HiCcompare` jointly normalizes two Hi-C datasets to remove biases between them. Then, it can detect significant differences between the datasets using a genomic distance-stratified permutation test. The novel concept of the MD plot, based on the commonly used MA plot or Bland-Altman plot, is the basis for these methods. The log **M**inus (difference) between chromatin interaction frequencies is plotted on the Y-axis. The genomic **D**istance is plotted on the X-axis. The MD plot allows for visualization, normalization, and comparing the differences between the Hi-C datasets in a distance-stratified manner.

The main functions are:

+ `hic_loess()` - performs joint `loess` normalization, minimizing global and local biases between Hi-C datasets
+ `hic_compare()` - performs the difference detection process to detect significant changes between Hi-C datasets and assist in comparative analysis

Several demo Hi-C datasets are also included in the package. Refer to the `HiCcompare` vignette for full usage instructions, `vignette("HiCcompare-vignette")`

## Installation

First, make sure you have all dependencies installed in R.

``` r
install.packages(c('dplyr', 'data.table', 'ggplot2', 'gridExtra', 
				   'mgcv', 'parallel', 'devtools'))

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(c("InteractionSet", "GenomicRanges", "IRanges", 
		   "BiocParallel", "QDNAseq", "GenomeInfoDbData"))			   
```

To install `HiCcompare` from Bioconductor, use the following commands.

``` r
# Bioconductor development version and GitHub Release contain major changes for difference detection
# it is recommended to use the GitHub release until the next Bioconductor update
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("HiCcompare")
library(HiCcompare)
```

Or, install the latest version of `HiCcompare` directly from the GitHub.

``` r
library(devtools)
install_github('dozmorovlab/HiCcompare', build_vignettes = TRUE)
library(HiCcompare)
```

## Usage

First, you will need to obtain some Hi-C data. Data is available from the sources listed in the overview, along with many others. You will need to extract the data and read it into R as either a 3 column sparse upper triangular matrix or a 7-column BEDPE file. For more details on data extraction, see the `HiCcompare`'s vignette.

Below is an example analysis using `HiCcompare`. The data in 3-column sparse upper triangular matrix format is loaded, and the first step is to create a `hic.table` object using the `create.hic.table()` function. Next, the two Hi-C matrices are jointly normalized using the `hic_loess()` function. Finally, difference detection can be performed using the `hic_compare()` function. The `hic_loess()` and `hic_compare()` functions will also produce an MD plot for visualizing the differences between the datasets. 

``` r
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

## Citations

Stansfield, John C., Kellen G. Cresswell, Vladimir I. Vladimirov, and Mikhail G. Dozmorov. [HiCcompare: An R-Package for Joint Normalization and Comparison of HI-C Datasets](https://doi.org/10.1186/s12859-018-2288-x). _BMC Bioinformatics_ 19, no. 1 (December 2018).

Stansfield, John C., Duc Tran, Tin Nguyen, and Mikhail G. Dozmorov. [R Tutorial: Detection of Differentially Interacting Chromatin Regions From Multiple Hi-C Datasets](https://doi.org/10.1002/cpbi.76). Current Protocols in Bioinformatics, May 2019

[HiCcompareWorkshop](https://github.com/mdozmorov/HiCcompareWorkshop) - "Detection of Differentially Interacting Chromatin Regions From Multiple Hi-C Datasets" workshop presented on Bioconductor 2020 conference

## Additional Vignettes

The `HiCcompare` paper included several supplemental files that showcase some of the methods' usage and reasoning. Below are the titles and brief descriptions of each of these vignettes along with links to the compiled `.pdf` and the source `.Rmd` files. 

[Normalization method comparison](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S1_File.pdf) - Comparison of several Hi-C normalization techniques to display the persistence of bias in individually normalized chromatin interaction matrices, and its effect on the detection of differential chromatin interactions. [Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S1_File.Rmd)

[Estimation of the IF power-law depencence](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S2_File.pdf) - Estimation of the power-law dependence between the $log_{10}-log_{10}$ interaction frequencies and distance between interacting regions. This vignette displays the reasoning behind using a power-law function to simulate the signal portion of Hi-C matrices. [Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S2_File.Rmd)

[Estimation of the SD power-law dependence](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S3_File.pdf) - Estimation of the power-law dependence between the $log_{10}-log_{10}$ SD of interaction frequencies and distance between interacting regions. This vignette displays the reasoning behind using a power-law function to simulate the noise component of Hi-C matrices. [Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S3_File.Rmd)

[Estimation of proportion of zeros](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S4_File.pdf) - Estimation of the dependence between the proportion of zeros and distance between interacting regions. This vignette shows the distribution of zeros in real Hi-C data. The results were used for modeling the proportion of zeros in simulated Hi-C matrices with a linear function. [Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S4_File.Rmd)

[Evaluation of difference detection in simulated data](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S5_File.pdf) - Extended evaluation of differential chromatin interaction detection analysis using simulated Hi-C data. Many different classifier performance measures are presented. Note: if trying to compile the source `.Rmd` this will take a long time to knit. [Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S5_File.Rmd)

[Evaluation of difference detection in real data](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S6_File.pdf) - Extended evaluation of differential chromatin interaction detection analysis using real Hi-C data. Many different classifier performance measures are presented. Note: if trying to compile the source `.Rmd` this will take a long time to knit. [Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S6_File.Rmd)

[loess at varying resolution](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S7_File.pdf) - Visualization of the `loess` joint normalization over varying resolutions. This vignette shows that increasing sparsity of Hi-C matrices with increasing resolution causes loess to become less useful for normalization at high resolutions. [Source](https://github.com/dozmorovlab/HiCcompare/raw/supplemental/supplemental_files/S7_File.Rmd)


## Contributions & Support

Suggestions for new features and bug reports are welcome. Please create a new [issue](https://github.com/dozmorovlab/HiCcompare/issues) for any of these or contact the author directly: [@jstansfield0](https://github.com/jstansfield0) (stansfieldjc@vcu.edu)

## Contributors

Authors: [@jstansfield0](https://github.com/jstansfield0) (stansfieldjc@vcu.edu) & [@mdozmorov](https://github.com/mdozmorov) (mikhail.dozmorov@vcuhealth.org)
