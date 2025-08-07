Wrapper for running CellChat
================
2025-08-07

## Introduction

CellChat is a popular R package for inferring cell-cell communication
from single-cell and/or spatial transcriptomics data. It offers
functions to visualize the analysis, and the github repository can be
found [**here**](https://github.com/jinworks/CellChat). CellChat has
interaction databases for three different species, of which we will use
HUMAN to find consensus with CellPhoneDB. Like CellPhoneDB, its database
contains only manually curated ligand-receptor interactions and the
subunit architecture of heteromeric complexes are considered. In
additon, CellChat considers co-factors that are involved in the
activation or inhibition of ligand-receptor interactions.

## Running CellChat directly from Seurat object

`run_cellchat()` takes as input a Seurat object with the normalized
counts in the data layer and annotated cell types stored in its
meta.data. For more details on `run_cellchat()`, users may refer
[**here**](../R/run_cellchat.R). CellChat offers different methods of
calculating average gene expression per cell type (default: triMean). As
described by the developers of CellChat, triMean is a robust mean method
given by \[(Q1 + 2M + Q3) / 4\] that returns fewer but stronger
interactions. If users would like to change the method for calculating
average expression, they may use the `type` parameter to manually
overwrite this.

`run_cellchat` takes in the same parameters as CellChat’s
`createCellChat()`, `subsetDB()`, `computeCommunProb()` and
`filterCommunication()`. Details on the respective functions can be
found by running help(<package_name>) in R.

## Tutorial

``` r
pacman::p_load_gh("nigelhojinker/CCCtools")
setwd(this.path::here())
rm(list = ls())

data("seu.NL")

seu.NL@meta.data %>% head()
```

    ##                     patient.id      labels
    ## S3_ATGAGGGAGTCTTGCA   Patient3   FBN1+ FIB
    ## S3_GTACGTACAAATTGCC   Patient3   FBN1+ FIB
    ## S3_CTCGTCAGTGTTGAGG   Patient3 Inflam. FIB
    ## S3_ATTCTACGTAATCGTC   Patient3   FBN1+ FIB
    ## S3_CTGCCTATCAATCACG   Patient3   FBN1+ FIB
    ## S3_TAGTGGTAGGATGCGT   Patient3   FBN1+ FIB

In `seu.Nl`, the metadata column for the annotated cell type labels is
“labels”. As such, we take “labels” as the input for the `group.by`
argument in `run_cellchat()`.

### Running CellChat on with all default parameters

This step typically takes a few minutes to run:

``` r
cellchat <- run_cellchat(seu.NL, group.by = "labels", assay = "RNA")

# [1] "Create a CellChat object from a Seurat object"
# The `meta.data` slot in the Seurat object is used as cell meta information 
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  APOE+ FIB, FBN1+ FIB, COL11A1+ FIB, Inflam. FIB, cDC1, cDC2, LC, Inflam. DC, TC, Inflam. TC, CD40LG+ TC, NKT 
# triMean is used for calculating the average gene expression per cell group. 
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
# The number of highly variable ligand-receptor pairs used for signaling inference is 440 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-08-07 15:54:08.823879]"
#   |=======================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-08-07 15:57:07.573151]"
# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  cDC1, Inflam. DC !   37.6% interactions are removed!
# CellChat analysis completed succesfully.
# Warning message:
# In createCellChat(obj, group.by = group.by, assay = assay) :
#   The 'meta' data does not have a column named `samples`. We now add this column and all cells are assumed to belong to `sample1`! 
```

### Running CellChat with filtered database

For users interested in only using a subset of CellChat’s database, they
may filter the database in accordance with the `subsetDB()` function by
CellChat. We provide an argument `subsetDB` that is set to FALSE by
default. To filter, simply set `subsetDB = TRUE` in `run_cellchat()` and
add in the corresponding arguments as per CellChat’s `subsetDB()`
function.

``` r
cellchat <- run_cellchat(seu.NL, group.by = "labels", assay = "RNA",
                         subsetDB = TRUE, search = c("Cell-Cell Contact"), key = c("annotation"))
```

For more than one level of filtering, `key` takes in a character vector
of the columns in CellChatDB and `search` has to take in a list
corresponding to the item(s) to filter in the keys.

``` r
cellchat <- run_cellchat(seu.NL, group.by = "labels", assay = "RNA",
                         subsetDB = TRUE, 
                         search = list(c("Cell-Cell Contact", "Secreted Signaling"), c("CellChatDB v1")),
                         key = c("annotation", "version"))
```

To remove all non-protein signaling interactions, since a majority of
non-protein signaling interactions consists of metabolic and synaptic
signaling:

``` r
cellchat <- run_cellchat(seu.NL, group.by = "labels", assay = "RNA",
                         subsetDB = TRUE)
```
