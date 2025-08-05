
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CCCtools

<!-- badges: start -->
<!-- badges: end -->

The goal of CCCtools is to provide functions to run CellChat v2 and
CellPhoneDB v5 as well as compare their outputs.

## Installation

You can install the development version of CCCtools:

``` r
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load_gh("Nigel-Ho23/CCCtools")
```

## Example data

`seu.NL` is the Seurat object contained 2,233 cells from non-lesional
human skin from four different patients. This dataset was created from
the data used originally in CellChat tutorial. See this [script for how
the data was generated](data-raw/demo_data.md).

``` r
library(CCCtools)

data(seu.NL)

seu.NL
#> Loading required package: SeuratObject
#> Warning: package 'SeuratObject' was built under R version 4.4.3
#> Loading required package: sp
#> Warning: package 'sp' was built under R version 4.4.2
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
#> An object of class Seurat 
#> 10353 features across 2233 samples within 1 assay 
#> Active assay: RNA (10353 features, 0 variable features)
#>  1 layer present: data
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
