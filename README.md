
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CCCtools

<!-- badges: start -->

<!-- badges: end -->

The goal of CCCtools is to provide functions to run [CellChat
v2](https://github.com/jinworks/CellChat) and [CellPhoneDB
v5](https://github.com/ventolab/CellphoneDB/tree/master) as well as
compare their outputs. Gokce was here!

## Installation

You can install the development version of CCCtools:

``` r
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load_gh("nigelhojinker/CCCtools")
```

## Example data

`seu.NL` is the Seurat object containing 2,233 cells from non-lesional
human skin from four different patients. This dataset was created from
the data used originally in CellChat tutorial. See this [script for how
the data was generated](data-raw/demo_data.md).

``` r
suppressPackageStartupMessages(library(CCCtools)) 
#> Warning: package 'purrr' was built under R version 4.5.1
#> Warning: package 'reticulate' was built under R version 4.5.1

data(seu.NL)

seu.NL
#> An object of class Seurat 
#> 10353 features across 2233 samples within 1 assay 
#> Active assay: RNA (10353 features, 0 variable features)
#>  1 layer present: data
```

## Tutorials

In **CCCtools**, we provide the functions `run_cellchat()` and
`run_cellphonedb()` for users to perform CellChat and CellPhoneDB
(method 2) analysis directly on their processed Seurat object.

Please refer to the following links to the documentations on running
these tools on your dataset:

- [Running CellChat on my Seurat object](../data-raw/Run_CellChat.md)
- [Running CellPhoneDB on my Seurat
  object](../data-raw/Run_CellPhoneDB.md)
