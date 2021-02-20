#' Preparation of a new huva_dataset with user-defined data
#'
#' @title generate_huva_dataset
#' @description This function provides a new dataset including expression data, sample annotation
#'              and possible metadata to run a huva experiment with novel data.
#' @param dataset_name character name for the new huva dataset.
#' @param data dataset input data.
#' @param annotation dataset annotation table.
#' @param metadata if NULL no metadata will be added.
#' @param combine LOGICAL, if TRUE the new dataset will be combined with the huva default dataset.
#' @return The sum of \code{x} and \code{y}.
#' @seealso run_huva_experiment
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' data_FG <- huva.db$FG500$data
#' anno_FG <- huva.db$FG500$anno
#' meta_FG <- huva.db$FG500$metadata
#'
#' new_huva.db <- generate_huva_dataset(dataset_name = "newdata",
#'                                      data = "data_FG",
#'                                      annotation = "anno_FG",
#'                                      metadata = "meta_FG",
#'                                      combine = FALSE)
#'
#'
#' @export
generate_huva_dataset <- function(dataset_name, data, annotation, metadata = NULL, combine = F) {

  container <- list()

  for (i in data) {
    container[[dataset_name]][["data"]] <- get(data)
  }

  for (j in annotation) {
    container[[dataset_name]][["anno"]] <- get(annotation)
  }

  for (k in metadata) {
    container[[dataset_name]][["metadata"]] <- get(metadata)
  }

  if (combine==T) {
    container <- c(huva_default_dataset, container)
  }

  class(container) <- "huva_dataset"

  return(container)
}

#' Plot heatmap object
#'
#' @title plot_HM
#' @description plot heatmap object
#' @param HM object containing output from the pheatmap function. Created by assigning output of
#'        pheatmap function to an object (i.e. pheatmap_object <- pheatmap(matrix))
#' @keywords plot_pheatmap
#' @return A heatmap created with pheatmap function
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' binned_dataset <- run_huva_experiment(data = huva.db,
#'                                      gene = "MYD88",
#'                                      quantiles = 0.10,
#'                                      gs_list = hallmarks_V7.2,
#'                                      summ = TRUE,
#'                                      datasets_list = NULL,
#'                                      adjust.method = "none")
#'
#' DE_huva <- get_DE_huva(huva_exp = binned_dataset,
#'                        study = "FG500",
#'                        dataset = "FG500_whole_blood")
#'
#' plot_HM(DE_huva$HM_FG500_whole_blood)
#'
#' @export
plot_HM <-function(HM){
  grid::grid.newpage()
  grid::grid.draw(HM$gtable)
}
