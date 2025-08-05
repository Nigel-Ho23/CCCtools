pacman::p_load(tidyverse, janitor, reticulate)
setwd(this.path::here())
rm(list = ls())

# Load your seurat object
seurat_obj <- readRDS("../data/tutorial_data.rds")


# Set up conda & python via reticulate ------------------------------------

# For users without a prior conda installation
# if(!file.exists(conda_binary())) install_miniconda()

# If conda installation is not found automatically:
options(reticulate.conda_binary = "C:/Users/hojkn/AppData/Local/miniforge3/condabin/conda.bat")

# Check that conda is available and detected
conda_binary()

# Create CellPhoneDB environment
# conda_create(envname = "cpdb", environment = "cpdb.yaml")

# Select python for analysis
python_binary <- conda_list() %>%
  filter(name == "cpdb") %>%
  pull(python) %>%
  normalizePath(winslash = "/")

python_binary

use_python(python_binary, required = TRUE)

pacman::p_load_gh("Nigel-Ho23/CCCtools")

test <- run_cellphonedb(seurat_obj,
                        condition = "condition",
                        toKeep    = "NL",
                        cpdbPath  = "../data/cellphonedb.zip",
                        metaPath  = "../data/NL_meta.tsv",
                        adataPath = "../data/NL.h5ad",
                        outPath   = "results/cpdb/NL")


help(run_cellphonedb)



