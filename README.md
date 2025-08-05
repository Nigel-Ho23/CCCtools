
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

`NL` is a Seurat object created from the example dataset in the the
CellChat tutorial. See the `create_NL.R` script in data-raw folder.

``` r
library(CCCtools)

hello("Nigel Ho")
#> [1] "Hello Nigel Ho"
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
