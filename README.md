# immApex
API for single-cell immune repertoire deep learning models

<img align="right" src="https://github.com/ncborcherding/immApex/blob/main/www/immApex_hex.png" width="305" height="352">

<!-- badges: start -->
[![R-CMD-check](https://github.com/ncborcherding/immApex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncborcherding/immApex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ncborcherding/immApex/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ncborcherding/immApex?branch=main)
<!-- badges: end -->

## Introduction

Single-cell sequencing is now a integral tool in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, 
like immune cell receptor profiling at the level of an individual cell. Towards this end, we developed the [scRepertoire](https://github.com/ncborcherding/scRepertoire) 
R package to assist in the interaction of immune receptor and gene expression sequencing. Further we developed models for embedding single-cell TCR sequences ([Trex](https://github.com/ncborcherding/Trex)) and BCR sequences ([Ibex](https://github.com/ncborcherding/Ibex)) using convolutional neural networks. **immApex** is the API for preparing the sequence data for the current and future models in the scRepertoire ecosystem. 

## System requirements 

immApex has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages - it is specifically designed to work with single-cell objects that have 
had BCR/TCRs added using [scRepertoire](https://github.com/ncborcherding/scRepertoire). immApex has been tested on OS X and Linux platforms.

**keras3** is necessary for immApex (this includes the set up of the tensorflow environment in R):

```r
##Install keras3
install.packages("keras3")

##Setting up Tensor Flow
library(reticulate)
use_condaenv(condaenv = "r-reticulate", required = TRUE)
library(keras3)
install_keras()
```

## Installation

To run immApex, open R and install immApex from github: 

```r
devtools::install_github("ncborcherding/immApex")
```

or via Bioconductor with the 3.20 release

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("immApex")
```

### IMGT usage

IMGT is used as a reference for gene names and sequence information can be accessed via ```getIMGT()```. Data from IMGT is under a CC BY-NC-ND 4.0 license. Please be aware that attribution is required for usage and it is the intent of IMGT to not allow derivative or commercial usage. 

## Usage/Demos

immApex should be able to be run in popular R-based single-cell workflows, including Seurat and Bioconductor/Single-Cell Experiment formats.

## Quick Start 

Check out this [vignette](https://www.borch.dev/uploads/screpertoire/articles/immapex) for a quick start tutorial. 

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/ncborcherding/immApex/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/ncborcherding/immApex/issues).

#### [Pull Requests](https://github.com/ncborcherding/immApex/pulls) are welcome for bug fixes, new features, or enhancements.
