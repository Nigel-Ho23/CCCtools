pacman::p_load(tidyverse, janitor, reticulate)
pacman::p_load_gh("Nigel-Ho23/CCCtools")
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


run_cellphonedb(obj = seurat_obj,
                group = "condition",
                toKeep    = "NL",
                cpdbPath  = "../data/cellphonedb.zip",
                metaPath  = "../data/NL_meta.tsv",
                adataPath = "../data/NL.h5ad",
                outPath   = "results/cpdb/NL")


help(run_cellphonedb)

ad <- import("anndata")
cpdb.analysis <- import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")

if (!is.null(group)) {
  to.keep = rownames(obj@meta.data)[obj@meta.data[[group]] == toKeep]
  obj <- subset(obj, cells = to.keep)
}

meta <- obj@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  select(barcode, labels) %>%
  write_delim(file = metaPath, delim = "\t") # save file for later use

normcounts <- obj[["RNA"]]$data

identical(meta$barcode, colnames(normcounts))

genes   <- data.frame(gene = rownames(normcounts)) %>%
  column_to_rownames("gene")

barcode <- data.frame(barcode = colnames(normcounts)) %>%
  column_to_rownames("barcode")

adata <- ad$AnnData(X   = normcounts,
                    obs = genes,
                    var = barcode)

adata <- adata$transpose()

adata$write_h5ad(adataPath)

dir.create(outPath, recursive = TRUE)

cpdb.analysis$call(
  cpdb_file_path = cpdbPath,
  meta_file_path = metaPath,
  counts_file_path = adataPath,
  counts_data = 'hgnc_symbol',
  # active_tfs_file_path = active_tf_path,
  # microenvs_file_path = microenvs_file_path,
  output_path = outPath
)



