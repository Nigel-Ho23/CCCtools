#' Run CellPhoneDB version 5 in R
#'
#' @description
#' Runs CellPhoneDBv5 on a Seurat object
#'
#' @param obj Seurat object with normalized counts in the data layer of the RNA assay
#' @param group Metadata column to filter by; this is for seurat objects containing cells from different disease status or time frames etc. (default = NULL)
#' @param toKeep Value within the metadata column of interest to retain for the analysis. This field is mandatory if group is not NULL. (default = NULL)
#' @param labels Metadata column name of cell type annotations
#' @param cpdbPath Path to 'cellphonedb.zip' file
#' @param metaPath Path to metadata file created internally - must end with .tsv
#' @param adataPath Path to AnnData object created internally - must end with .h5ad
#' @param outPath Path to store all output .txt files from CellPhoneDB analysis
#' @param counts_data see Value below for link to description by CellPhoneDB
#' @param active_tfs_file_path see Value below for link to description by CellPhoneDB
#' @param microenvs_file_path see Value below for link to description by CellPhoneDB
#' @param score_interactions see Value below for link to description by CellPhoneDB
#' @param iterations see Value below for link to description by CellPhoneDB
#' @param threshold see Value below for link to description by CellPhoneDB
#' @param threads see Value below for link to description by CellPhoneDB
#' @param debug_seed see Value below for link to description by CellPhoneDB
#' @param result_precision see Value below for link to description by CellPhoneDB
#' @param pvalue see Value below for link to description by CellPhoneDB
#' @param subsampling see Value below for link to description by CellPhoneDB
#' @param subsampling_log see Value below for link to description by CellPhoneDB
#' @param subsampling_num_pc see Value below for link to description by CellPhoneDB
#' @param subsampling_num_cells see Value below for link to description by CellPhoneDB
#' @param separator see Value below for link to description by CellPhoneDB
#' @param debug see Value below for link to description by CellPhoneDB
#' @param output_suffix see Value below for link to description by CellPhoneDB
#'
#' @returns Folder path to CellPhoneDB analysis outputs
#'
#' @returns For arguments from counts_data to output_suffix, please refer to https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T1_Method2.ipynb
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Runs CellPhoneDB method 2 with all default parameters
#' run_cellphonedb(
#'   obj = seurat_object,
#'   group  = "condition",           # if seurat object contains multiple conditions to filter, input column name in metadata
#'   toKeep = "treated",             # mandatory if group is not NULL; input condition to be filtered
#'   labels = "celltype",
#'   cpdbPath  = "data/cellphonedb.zip",
#'   metaPath  = "data/meta.tsv",    # internally created metadata file will be stored here
#'   adataPath = "data/data.h5ad",   # internally created AnnData object will be stored here
#'   outPath   = "output/results")
#'
#' ## For users interested in customising their CellPhoneDB run
#' run_cellphonedb(seu.NL,
#'   labels = "labels",
#'   cpdbPath = "../data/cellphonedb.zip",
#'   metaPath = "../test/testagain/NL_meta.tsv",
#'   adataPath = "../test/testagain/NL_normcounts.h5ad",
#'   outPath = "results/cpdb/NL",
#'   iterations = 100,
#'   threshold = 0.2,
#'   threads = 5,
#'   debug_seed = 42,
#'   result_precision = 5,
#'   score_interactions = TRUE)
#' }

run_cellphonedb <- function (obj, group.by = NULL, toKeep = NULL, labels, cpdbPath,
                             metaPath, adataPath, counts_data = "hgnc_symbol", active_tfs_file_path = NULL, microenvs_file_path = NULL,
                             score_interactions = FALSE, iterations = 1000, threshold = 0.1, threads = 4, debug_seed = -1,
                             result_precision = 3, pvalue = 0.05, subsampling = FALSE, subsampling_log = FALSE, subsampling_num_pc = 100,
                             subsampling_num_cells = 1000, separator = "|", debug = FALSE, outPath, output_suffix = NULL) {

  ## Convert arguments to integers for CellPhoneDB python function
  iterations <- as.integer(iterations)
  threads <- as.integer(threads)
  debug_seed <- as.integer(debug_seed)
  result_precision <-as.integer(result_precision)
  subsampling_num_pc <- as.integer(subsampling_num_pc)
  subsampling_num_cells <- as.integer(subsampling_num_cells)

  ## Load required modules
  ad <- import("anndata")
  cpdb.analysis <- import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")

  ## Input file creation
  if (!is.null(group.by)) {
    to.keep <- rownames(obj@meta.data)[obj@meta.data[[group.by]] == toKeep]
    obj <- subset(obj, cells = to.keep)
  }

  for (path in c(metaPath, adataPath)) {
    if (!dir.exists(dirname(path))) {
      dir.create(dirname(path), recursive = TRUE)
    }
  }

  meta <- obj@meta.data %>% rownames_to_column(var = "barcode") %>%
    select(barcode, !!rlang::sym(labels)) %>% write_delim(file = metaPath, delim = "\t")

  normcounts <- obj[["RNA"]]$data

  identical(meta$barcode, colnames(normcounts))

  genes <- data.frame(gene = rownames(normcounts)) %>%
    column_to_rownames("gene")

  barcode <- data.frame(barcode = colnames(normcounts)) %>%
    column_to_rownames("barcode")

  ## Create AnnData object
  adata <- ad$AnnData(X = normcounts, obs = genes, var = barcode)
  adata <- adata$transpose()
  adata$write_h5ad(adataPath)

  if (!dir.exists(outPath)) {
    dir.create(outPath, recursive = TRUE)
  }

  ## Run analysis
  cpdb.analysis$call(cpdb_file_path = cpdbPath,
                     meta_file_path = metaPath,
                     counts_file_path = adataPath,
                     counts_data = counts_data,
                     active_tfs_file_path = active_tfs_file_path,
                     microenvs_file_path = microenvs_file_path,
                     score_interactions = score_interactions,
                     iterations = iterations,
                     threshold = threshold,
                     threads = threads,
                     debug_seed = debug_seed,
                     result_precision = result_precision,
                     pvalue = pvalue,
                     subsampling = subsampling,
                     subsampling_log = subsampling_log,
                     subsampling_num_pc = subsampling_num_pc,
                     subsampling_num_cells = subsampling_num_cells,
                     separator = separator,
                     debug = debug,
                     output_path = outPath,
                     output_suffix = output_suffix)

  message("CellPhoneDB analysis completed successfully.")

}
