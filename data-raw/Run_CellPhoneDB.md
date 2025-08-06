Running CellPhoneDB in R
================
2025-08-06

# Introduction

CellPhoneDB is a popular **Python** package developed by [Teichmann
Lab](https://www.teichlab.org/), and currently managed and further
developed by the [Vento-Tormo Lab](https://ventolab.org/) from
CellPhoneDB version 3. It is a useful cell-cell communication tool to
directly look up interactions of interest from single-cell multiomics
data, with their database of only manually curated and reviewed HUMAN
interactions.

A key feature of CellPhoneDB is its inherent function of accounting for
heteromeric complexes by decomplexifying to their subunits and ensuring
that expression of all subunits have to be non-zero and above a
user-defined threshold to be considered present in the interaction. This
accurately recapitulates how interactions can occur in vivo. As such, we
thought to incorporate CellPhoneDB with another popular cell-cell
communication tool, CellChat, that has the same feature.

# Setting up required environment

As mentioned in the **Introduction**, CellPhoneDB is a Python package
and currently, does not have an R equivalent. Here, we provide a wrapper
function `run_cellphonedb()` for running CellPhoneDB completely within
R. To do this, we require the use of the `reticulate` package (a
dependency in CCCtools) and provide users with the `cpdb.yaml` file for
creating the CellPhoneDB python environment in R. In addition to the
`cpdb.yaml` file, the CellPhoneDB interaction database `cellphonedb.zip`
is also **required** for running CellPhoneDB analysis; you may find the
environment file and database in the same directory for download
[**here**](../data/).

## CellPhoneDB conda environment in R

### R set-up

``` r
pacman::p_load_gh("nigelhojinker/CCCtools")
setwd(this.path::here())
rm(list = ls())
```

### Set up conda & python via reticulate

For users **without** a prior conda installation, you may run the line
below for installation of miniconda:

``` r
if(!file.exists(conda_binary())) install_miniconda()
```

For users **with** conda installed, we recommend running the line below
to ensure that the correct conda is selected and used by `reticulate`.
Of course, users should input the path to conda in their respective
local machines.

``` r
options(reticulate.conda_binary = "C:/Users/hojkn/AppData/Local/miniforge3/condabin/conda.bat")
```

Once installation and conda selection is complete, we can confirm this
by running:

``` r
conda_binary()
```

    ## [1] "C:/Users/hojkn/AppData/Local/miniforge3/condabin/conda.bat"

Next, we have to create the CellPhoneDB conda environment with the
necessary Python modules for running the analysis. This is done **only
once:**

``` r
conda_create(envname = "cpdb", environment = "/path/to/cpdb.yaml")
```

Once the cpdb environment has been created, we have to select its python
interpreter:

``` r
python_binary <- conda_list() %>%
  filter(name == "cpdb") %>%
  pull(python) %>%
  normalizePath(winslash = "/")

python_binary
```

    ## [1] "C:/Users/hojkn/AppData/Local/miniforge3/envs/cpdb/python.exe"

``` r
use_python(python_binary, required = TRUE)
```

## Running CellPhoneDB in R

Now, we are ready to run CellPhoneDB in R.

We provide the `run_cellphonedb()` function that can give us output from
CellPhoneDB analysis (method 2). This function takes in a Seurat object
with normalized counts in the data layer of the RNA assay, and annotated
cell type labels in the meta.data. Here is a breakdown of the other
parameters in `run_cellphonedb()`:

- labels: meta.data column name for the annotated cell type labels
- cpdbPath: /path/to/cellphonedb.zip
- metaPath: /path/to/meta.tsv; path to the tsv file and the tsv file
  itself will be created internally, there is no need to make the file
  path before running this function
- adataPath /path/to/normcounts.h5ad; function creates an AnnData object
  internally and save it to the user defined path. Similar to metaPath,
  users do not have to create the file path before running the function
- outPath: /path/to/outputs; five .txt files will be generated in every
  run and saved in this directory (creation prior to run is not needed)

``` r
# Load your Seurat object
data(seu.NL)

seu.NL

seu.NL@meta.data %>% head()

run_cellphonedb(obj = seu.NL,
                labels = "labels",
                cpdbPath  = "../data/cellphonedb.zip",
                metaPath  = "../data/NL_meta.tsv",
                adataPath = "../data/NL.h5ad",
                outPath   = "results/test/NL")
```
