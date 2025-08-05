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
conda_create(envname = "cpdb", environment = "cpdb.yaml")

# Select python for analysis
python_binary <- conda_list() %>%
  filter(name == "cpdb") %>%
  pull(python) %>%
  normalizePath(winslash = "/")

python_binary

use_python(python_binary, required = TRUE)


# Run CellPhoneDB ---------------------------------------------------------

ad <- import("anndata")
cpdb.analysis <- import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")

NL.meta <- seurat_obj@meta.data %>%
  filter(condition == "NL") %>%
  rownames_to_column(var = "barcode") %>%
  select(barcode, labels) %>%
  write_delim(file = "../data/NL_meta.tsv", delim = "\t") # save file for later use

NL <- subset(seurat_obj, condition == "NL")

NL.normcounts <- NL[["RNA"]]$data

identical(NL.meta$barcode, colnames(NL.normcounts))

NL.genes   <- data.frame(gene = rownames(NL.normcounts)) %>%
  column_to_rownames("gene")

NL.barcode <- data.frame(barcode = colnames(NL.normcounts)) %>%
  column_to_rownames("barcode")

adata <- ad$AnnData(X   = NL.normcounts,
                    obs = NL.genes,
                    var = NL.barcode)

adata <- adata$transpose()

adata$write_h5ad("../data/NL.h5ad")

# Main database containing the interaction
download.file("https://raw.github.com/ventolab/cellphonedb-data/refs/heads/master/cellphonedb.zip",  "../data/cellphonedb.zip")

cpdb_file_path   <- "../data/cellphonedb.zip"
meta_file_path   <- "../data/NL_meta.tsv"
counts_file_path <- "../data/NL.h5ad"
out_path         <- "results/cpdb/NL"

dir.create(out_path, recursive = TRUE)

cpdb_results <- cpdb.analysis$call(
  cpdb_file_path = cpdb_file_path,
  meta_file_path = meta_file_path,
  counts_file_path = counts_file_path,
  counts_data = 'hgnc_symbol',
  # active_tfs_file_path = active_tf_path,
  # microenvs_file_path = microenvs_file_path,
  output_path = out_path
)

"NL"
"../data/cellphonedb.zip"
"../data/NL_meta.tsv"
"../data/NL.h5ad"
"results/cpdb/NL"
