#' Demo dataset prepared from CellChat tutorial data
#'
#' The dataset used in the main CellChat tutorial
#' after quality control filtering (min.cells = 75, min.features = 200)
#'
#' @format A Seurat object with 10,353 genes and 2,233 cells. The meta data slot includes:
#' \describe{
#'   \item{patient.id}{Patient identifier}
#'   \item{labels}{Cell type annotation provided by the CellChat authors}
#' }
#'
#' @examples
#'
#' data(seu.NL)
#'
#' seu.NL@meta.data %>%
#'   tabyl(labels)
#'
#' @source
#' Data was downloaded from https://doi.org/10.6084/m9.figshare.13520015.v1 and processed further (see https://github.com/Nigel-Ho23/CCCtools/tree/main/data-raw/demo_data.md for details).
"seu.NL"
