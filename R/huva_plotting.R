#' Plot selected genes from an huva_experiment class object
#'
#' @title plot_binned_gene
#' @description This function plots the expression of selected genes in a huva experiment low- and high- groups across
#'              the different provided datasets.
#' @param goi name of the gene/s of interest
#' @param huva_experiment huva_experiment class object
#' @return list of ggplot objects, each plot recapitulates gene expression of selected gene/s in the present datasets.
#' @import ggplot2
#' @import reshape2
#' @seealso run_huva_experiment
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
#' plot_binned <- plot_binned_gene(goi = c("CRELD1", "CRELD2", "TCF7", "RCAN3"),
#'                                 huva_experiment = binned_dataset)
#'
#' @export
plot_binned_gene <- function(goi, huva_experiment) {

  if (class(huva_experiment)!= "huva_experiment") {
    print("use huva_experiment class object for reliable results")
  }

  container <- list()

  for (i in names(huva_experiment)) {

    for (j in names(huva_experiment[[i]][["data"]])) {

      if (sum(goi %in% rownames(huva_experiment[[i]][["data"]][[j]]))>=1) {

        tmp <- as.data.frame(huva_experiment[[i]][["data"]][[j]])
        tmp <- t(tmp[goi,])

        #q <- paste("anno", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")

        tmp <- merge(huva_experiment[[i]][["anno"]][[j]][, c("Row.names", "group")], tmp, by.x = "Row.names", by.y = "row.names")
        tmp <- suppressMessages(melt(tmp, value.name = "expression", variable.name = "gene"))

        tmp <- ggplot(tmp, aes(x = gene, y = expression, fill = group))+
          geom_boxplot() + xlab("") + ylab("expression") + ggtitle(paste(unlist(strsplit(j, "_"))[-1], collapse = "_")) +
          theme_minimal() + theme(aspect.ratio = 1/length(goi)) + stat_compare_means(method = "t.test") +
          scale_fill_manual(values = c("#c93918", "#0e2e99"))

        container[[j]] <- tmp

      }

    }

  }

  return(container)

}
