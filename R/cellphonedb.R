#' Run CellPhoneDB version 5 in R
#'
#' @param obj Seurat object
#' @param ...
#'
#' @returns Folder path
#' @export
#'
run_cellphonedb <- function(obj, condition = NULL, toKeep = NULL, cpdbPath, metaPath, adataPath, outPath) {
  ad <- import("anndata")
  cpdb.analysis <- import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")

  if (is.null(condition)) {
    obj <- obj
  }
  else {
    obj <- subset(obj, obj@meta.data[[condition]] == toKeep)
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

  cpdb_results <- cpdb.analysis$call(
    cpdb_file_path = cpdbPath,
    meta_file_path = metaPath,
    counts_file_path = adataPath,
    counts_data = 'hgnc_symbol',
    # active_tfs_file_path = active_tf_path,
    # microenvs_file_path = microenvs_file_path,
    output_path = outPath
  )

  print("done")

}
