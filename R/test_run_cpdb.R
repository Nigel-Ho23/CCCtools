pacman::p_load_gh("nigelhojinker/CCCtools")
setwd(this.path::here())
rm(list = ls())

# Load your seurat object
data(seu.NL)


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


# Run CellPhoneDB in R ----------------------------------------------------

help(run_cellphonedb)

run_cellphonedb(obj = seu.NL,
                cpdbPath  = "../data/cellphonedb.zip",
                metaPath  = "../test/testagain/NL_meta.tsv",
                adataPath = "../test/testagain/NL.h5ad",
                outPath   = "results/test/testing/NL")
