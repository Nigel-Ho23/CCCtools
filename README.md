
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

`seu.NL` is the Seurat object contained 2,233 cells from non-lesional
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
