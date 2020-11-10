#' Preparation of a new huva_dataset with independent data
#'
#' @param dataset_name charachter name for the new huva dataset generated
#' @param data dataset data
#' @param annotation dataset annotation
#' @param metadata if NULL no metadata will be added
#' @param combine LOGICAL, if TRUE the new dataset will be combined with the huva default dataset
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' # This will generate a combined dataset with new data and the data provided in the huva package
#'
#' 2+2
#'
#' @export
generate_huva_dataset <- function(dataset_name, data, annotation, metadata=NULL, combine=F) {

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
#' @param HM object containing output from the pheatmap function. Created by assigning output of
#' pheatmap function to an object (i.e. pheatmap_object <- pheatmap(matrix))
#' @keywords plot_pheatmap
#' @export
#' @return A heatmap created by pheatmap
#' @examples
#' 2+2

plot_HM <-function(HM){
  grid::grid.newpage()
  grid::grid.draw(HM$gtable)
}
