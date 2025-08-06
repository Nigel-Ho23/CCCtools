
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
library(CCCtools)
#> Loading required package: Seurat
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
#> Loading required package: tidyverse
#> Warning: package 'purrr' was built under R version 4.5.1
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.2     ✔ tibble    3.3.0
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
#> ✔ purrr     1.1.0
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
#> Loading required package: janitor
#> 
#> 
#> Attaching package: 'janitor'
#> 
#> 
#> The following objects are masked from 'package:stats':
#> 
#>     chisq.test, fisher.test
#> 
#> 
#> Loading required package: CellChat
#> 
#> Loading required package: igraph
#> 
#> 
#> Attaching package: 'igraph'
#> 
#> 
#> The following objects are masked from 'package:lubridate':
#> 
#>     %--%, union
#> 
#> 
#> The following objects are masked from 'package:dplyr':
#> 
#>     as_data_frame, groups, union
#> 
#> 
#> The following objects are masked from 'package:purrr':
#> 
#>     compose, simplify
#> 
#> 
#> The following object is masked from 'package:tidyr':
#> 
#>     crossing
#> 
#> 
#> The following object is masked from 'package:tibble':
#> 
#>     as_data_frame
#> 
#> 
#> The following object is masked from 'package:Seurat':
#> 
#>     components
#> 
#> 
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> 
#> 
#> The following object is masked from 'package:base':
#> 
#>     union
#> 
#> 
#> Loading required package: reticulate
#> Warning: package 'reticulate' was built under R version 4.5.1

data(seu.NL)

seu.NL
#> An object of class Seurat 
#> 10353 features across 2233 samples within 1 assay 
#> Active assay: RNA (10353 features, 0 variable features)
#>  1 layer present: data
```
