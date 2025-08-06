#' Run CellPhoneDB version 5 in R
#'
#' @param obj Seurat object with normalized counts in the data layer of the RNA assay
#' @param group Metadata column to filter by; this is for seurat objects containing cells from different disease status or time frames etc. (default = NULL)
#' @param toKeep Value within the metadata column of interest to retain for the analysis. This field is mandatory if group is not NULL. (default = NULL)
#' @param cpdbPath Path to 'cellphonedb.zip' file
#' @param metaPath Path to metadata file created internally - must end with .tsv
#' @param adataPath Path to AnnData object created internally - must end with .h5ad
#' @param outPath Path to store all output .txt files from CellPhoneDB analysis
#'
#' @returns Folder path to CellPhoneDB analysis outputs
#' @export
#'
#' @examples
#' \dontrun{
#' run_cellphonedb(
#'   obj = seurat_object,
#'   group  = "condition",           # if seurat object contains multiple conditions to filter, input column name in metadata
#'   toKeep = "treated",             # mandatory if group is not NULL; input condition to be filtered
#'   cpdbPath  = "data/cellphonedb.zip",
#'   metaPath  = "data/meta.tsv",    # internally created metadata file will be stored here
#'   adataPath = "data/data.h5ad",   # internally created AnnData object will be stored here
#'   outPath   = "output/results"
#' )
#' }

run_cellphonedb <- function(obj, group = NULL, toKeep = NULL, cpdbPath, metaPath, adataPath, outPath) {
  ad <- import("anndata")
  cpdb.analysis <- import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")

  if (!is.null(group)) {
    to.keep <- rownames(obj@meta.data)[obj@meta.data[[group]] == toKeep]
    obj <- subset(obj, cells = to.keep)
  }

  meta <- obj@meta.data %>%
    rownames_to_column(var = "barcode") %>%
    select(barcode, labels) %>%
    write_delim(file = metaPath, delim = "\t")

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

  if (!dir.exists(outPath)) {
    dir.create(outPath, recursive = TRUE)
    }

  cpdb.analysis$call(
    cpdb_file_path = cpdbPath,
    meta_file_path = metaPath,
    counts_file_path = adataPath,
    counts_data = 'hgnc_symbol',
    # active_tfs_file_path = active_tf_path,
    # microenvs_file_path = microenvs_file_path,
    output_path = outPath
  )

  message("CellPhoneDB analysis completed successfully.")

}
