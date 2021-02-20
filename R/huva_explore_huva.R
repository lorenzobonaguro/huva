#' get_expr_huva
#'
#' @title get_expr_huva
#' @description The function get_expr_huva provides the expression table resulting from reference huva experiment data.frame.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @return expression table of binned samples as data frame object.
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
#' expr_huva <- get_expr_huva(huva_exp = binned_dataset,
#'                            study = "FG500",
#'                            dataset = "FG500_whole_blood")
#'
#' @export
get_expr_huva <- function(huva_exp, study, dataset = NULL) {

  if (class(huva_exp)!= "huva_experiment") {
    error("Invalid class of huva experiment used as input, it can lead to incorrect results.")
  }

  if (! study %in% names(huva_exp)) {
    error("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["data"]])) {
      error("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    container <- huva_exp[[study]][["data"]][[dataset]]
  }
  else {
    container <- huva_exp[[study]][["data"]]
  }

  return(container)
}

#' get_anno_huva
#'
#' @title get_anno_huva
#' @description The function get_anno_huva produces a filtered annotation table including only donors belonging to the specified huva experiment groups.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @return annotation table of binned samples as data frame object.
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
#' anno_huva <- get_anno_huva(huva_exp = binned_dataset,
#'                            study = "FG500")
#'
#' @export
get_anno_huva <- function(huva_exp, study, dataset = NULL) {

  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["anno"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    container <- huva_exp[[study]][["anno"]][[dataset]]
  }
  else {
    container <- huva_exp[[study]][["anno"]]
  }

  return(container)
}

#' get_meta_huva
#'
#' @title get_meta_huva
#' @description The function get_meta_huva returns the metadata table of the specified huva_experiment groups.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @return annotation table of binned samples.
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
#' meta_huva <- get_meta_huva(huva_exp = binned_dataset,
#'                            study = "FG500")
#'
#' @export
get_meta_huva <- function(huva_exp, study, dataset = NULL) {

  if (class(huva_exp)!= "huva_experiment") {
    error("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    error("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (! "metadata" %in% names(huva_exp[[study]])) {
    error("No metadata available for this study.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["metadata"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    container <- huva_exp[[study]][["metadata"]][[dataset]]
  }
  else {
    container <- huva_exp[[study]][["metadata"]]
  }

  return(container)
}

#' get_DE_huva
#'
#' @title get_DE_huva
#' @description The function get_DE_huva returns the list of differentially expressed genes between groups of the specified huva experiment.In this function, the p-value (pval) and logFC (logFC) cutoffs can also be specified (default is logFC=1). Along with providing the DE genes table, this function also performs the Principal Component Analysis (PCA), which can be plotted over the desidered components (PC, default = c("PC1", "PC2")).
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @param cluster_col heatmap column clustering option, default is TRUE.
#' @param pval p-value cutoff used for the definition of differential expression.
#' @param logFC log2 fold change cutoff for the definition of differential expression.
#' @param PC character vector defining the principal components to be plotted, default is c("PC1", "PC2").
#' @return table of differentially expressed genes in the comparison between binned samples.
#' @seealso run_huva_experiment
#' @import pheatmap
#' @import grDevices
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
#' #Plot PCA
#' DE_huva$PCA_FG500_PBMC
#'
#' #Plot Histogram
#' DE_huva$plot_FG500_PBMC
#'
#' #Plot HM of DE genes
#' plot_HM(DE_huva$HM_FG500_PBMC)
#'
#' @export
get_DE_huva <- function(huva_exp, study, dataset = NULL, pval = 0.05, logFC = 0, cluster_col = T, PC = c("PC1", "PC2")) {
  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }
  if (! study %in% names(huva_exp)) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }
  if (!is.null(dataset)) {
    if (! dataset %in% names(huva_exp[[study]][["DE_genes"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }
    container <- list()
    data <- huva_exp[[study]][["DE_genes"]][[dataset]]
    data <- data[data$adj.P.Val < pval, ]
    data <- data[abs(data$logFC) > logFC, ]
    data$direction <- ifelse(data$logFC < 0, "down", "up")
    container[[dataset]] <- data
    container[[paste("plot", dataset, sep = "_")]] <- ggplot(as.data.frame(table(data$direction)), aes(x = Var1, y = Freq, fill = Var1))+
      geom_bar(stat = "identity", colour="black", size=0.5) + theme_minimal() +
      geom_text(stat = "identity", aes(label=Freq), vjust= -0.2) +
      theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_text(size = 12), axis.title = element_text(size = 14)) +
      xlab("") + ylab("# of DE genes") + ggtitle(paste("Number of ", dataset, sep = "")) + scale_fill_manual(values = c("#32a9d1", "#c92a2a"))
    # HM DE_genese
    hm_data <- huva_exp[[study]][["data"]][[dataset]][rownames(data),]
    hm_meta <- huva_exp[[study]][["anno"]][[dataset]]
    hm_meta <- hm_meta[order(hm_meta$group),]
    hm_data <- hm_data[,hm_meta$Row.names]
    rownames(hm_meta) <- hm_meta$Row.names
    container[[paste("HM", dataset, sep = "_")]] <- pheatmap(mat = hm_data,
                                                             scale = "row",
                                                             cluster_rows = T,
                                                             cluster_cols = cluster_col,
                                                             annotation = hm_meta[, c(2,3)],
                                                             annotation_colors = list(group = c(high = "#c93918", low = "#0e2e99"),
                                                                                      expression = heat.colors(20, alpha = 1)), silent = T)
    # PCA
    pca_data <- huva_exp[[study]][["data"]][[dataset]]
    pca_meta <- hm_meta
    pc <- prcomp(t(pca_data))$x
    pc <- merge(pc, pca_meta, by = "row.names")
    container[[paste("PCA", dataset, sep = "_")]] <- ggplot(pc, aes(x =.data[[PC[1]]], y =.data[[PC[2]]], colour = group))+
      geom_point(size = 5) +
      scale_colour_manual(values = c("#c93918", "#0e2e99")) +
      theme_bw() + theme(aspect.ratio = 1) +
      ggtitle(paste("PCA", dataset, sep = "_"))
  }
  else {
    container <- huva_exp[[study]]
    container_tmp <- list()
    for (i in names(container[["DE_genes"]])) {
      print(i)
      data <- container[["DE_genes"]][[i]]
      data <- data[data$adj.P.Val < pval, ]
      data <- data[abs(data$logFC) > logFC, ]
      data$direction <- ifelse(data$logFC < 0, "down", "up")
      container_tmp[[i]] <- data
      container_tmp[[paste("plot", i, sep = "_")]] <- ggplot(as.data.frame(table(data$direction)), aes(x = Var1, y = Freq, fill = Var1))+
        geom_bar(stat = "identity", colour = "black", size = 0.5) + theme_minimal() +
        geom_text(stat = "identity", aes(label = Freq), vjust = -0.2) +
        theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_text(size = 12), axis.title = element_text(size = 14)) +
        xlab("") + ylab("# of DE genes") + ggtitle(paste("Number of ", i, sep = "")) + scale_fill_manual(values = c("#32a9d1", "#c92a2a"))
      hm_data <- container[["data"]][[i]][rownames(data),]
      hm_meta <- container[["anno"]][[i]]
      hm_meta <- hm_meta[order(hm_meta$group),]
      hm_data <- hm_data[,hm_meta$Row.names]
      rownames(hm_meta) <- hm_meta$Row.names
      container_tmp[[paste("HM", i, sep = "_")]] <- pheatmap(mat = hm_data,
                                                             scale = "row",
                                                             cluster_rows = T,
                                                             cluster_cols = cluster_col,
                                                             annotation = hm_meta[, c(2,3)],
                                                             annotation_colors = list(group = c(high = "#c93918", low = "#0e2e99"),
                                                                                      expression = heat.colors(20, alpha = 1)), silent = T)
      pca_data <- container[["data"]][[i]]
      pca_meta <- hm_meta
      pc <- prcomp(t(pca_data))$x
      pc <- merge(pc, pca_meta, by="row.names")
      container_tmp[[paste("PCA", i, sep = "_")]] <- ggplot(pc, aes( x =.data[[PC[1]]], y =.data[[PC[2]]], colour=group))+
        geom_point(size = 5) +
        scale_colour_manual(values = c("#c93918", "#0e2e99")) +
        theme_bw() + theme(aspect.ratio = 1) +
        ggtitle(paste("PCA", i, sep = "_"))
    }
    container <- container_tmp
  }
  return(container)
}

#' get_rank_huva
#'
#' @title get_rank_huva
#' @description The function get_rank_huva extracts the log2FC-ranked gene list from the comparison between the "low" and "high" groups of the specified huva experiment.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @param n_top_genes number of top down- or up- regulated genes to be displayed.
#' @return log2FC-ranked gene list.
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
#' rank_huva <- get_rank_huva(huva_exp = binned_dataset,
#'                            study = "ImmVar",
#'                            dataset = NULL,
#'                            n_top_genes = 5)
#'
#' @export
get_rank_huva <- function(huva_exp, study, dataset = NULL, n_top_genes = 20) {

  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["Rank_genelist"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    rank <- huva_exp[[study]][["Rank_genelist"]][[dataset]]
    rank_top <- rbind(data.frame(delta_exp = rank[1:n_top_genes], gene_name = names(rank[1:n_top_genes]), group = "UP"),
                      data.frame(delta_exp = rank[(length(rank)-n_top_genes+1):length(rank)], gene_name = names(rank[(length(rank)-n_top_genes+1):length(rank)]), group = "DOWN"))

    rank_top$gene_name <- factor(rank_top$gene_name, levels = rank_top$gene_name)

    rank_top$group <- factor(rank_top$group, levels = c("UP", "DOWN")) # Bug fix of version 0.1.5

    plot <- ggplot(rank_top, aes(x = gene_name, y = delta_exp, fill = group)) +
      geom_point(size = 6, shape = 21) + facet_wrap(~group, scales = "free") + coord_flip() +
      theme_minimal() + theme(aspect.ratio = 2) + xlab("gene name") + ylab("delta expression low vs high") + scale_fill_manual(values = c("#c92a2a", "#32a9d1"))

    container <- list()

    container[[dataset]] <- rank_top

    container[[paste("plot", dataset, sep = "_")]] <- plot

  }

  else {

    container <- huva_exp[[study]][["Rank_genelist"]]

    container_tmp <- list()

    for (i in names(container)) {

      rank <- container[[i]]
      rank_top <- rbind(data.frame(delta_exp = rank[1:n_top_genes], gene_name = names(rank[1:n_top_genes]), group = "UP"),
                        data.frame(delta_exp = rank[(length(rank)-n_top_genes+1):length(rank)], gene_name = names(rank[(length(rank)-n_top_genes+1):length(rank)]), group = "DOWN"))

      rank_top$gene_name <- factor(rank_top$gene_name, levels = rank_top$gene_name)

      rank_top$group <- factor(rank_top$group, levels = c("UP", "DOWN")) # Bug fix of version 0.1.5

      plot <- ggplot(rank_top, aes(x = gene_name, y = delta_exp, fill = group)) +
        geom_point(size = 6, shape = 21) + facet_wrap(~group, scales = "free") + coord_flip() +
        theme_minimal() + theme(aspect.ratio = 2) + xlab("gene name") + ylab("delta expression low vs high") + scale_fill_manual(values = c("#c92a2a", "#32a9d1"))

      container_tmp[[i]] <- rank_top

      container_tmp[[paste("plot", i, sep = "_")]] <- plot

    }

    container <- container_tmp

  }

  return(container)

}

#' get_gsea_huva
#'
#' @title get_gsea_huva
#' @description Gene set enrichment analysis (GSEA) is performed on user-defined huva experiment ranked lists. By default, the huva experiment works on hallmark gene sets.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @param int_plot default is TRUE and provides an interactive GSEA Volcano plot.
#' @return return the GSEA table as a data frame and a Volcano plot (NES in x-axis, -log10pval in y-axis).
#' @seealso run_huva_experiment
#' @importFrom plotly plot_ly
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
#' gsea_huva <- get_gsea_huva(huva_exp = binned_dataset,
#'                            study = "FG500")
#'
#' gsea_huva$plot_FG500_whole_blood
#'
#' gsea_huva$int_plot_FG500_whole_blood
#'
#' @export
get_gsea_huva <- function(huva_exp, study, dataset = NULL, int_plot = T){

  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["gsea"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    gsea <- huva_exp[[study]][["gsea"]][[dataset]]

    plot <- ggplot(gsea, aes(x = NES, y = -log10(pval), size = -log10(pval), alpha = -log10(pval)))+
      geom_point() + theme_bw() + theme(aspect.ratio = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_x_continuous(limits = c(-max(abs(gsea$NES)), max(abs(gsea$NES))))

    container <- list()

    container[[dataset]] <- gsea
    container[[paste("plot", dataset, sep = "_")]] <- plot

    if (int_plot==T) {
      i_plot <- plot_ly(type = "scatter",
                        mode = 'markers',
                        data = gsea,
                        x = ~NES, y = ~(-log10(pval)),
                        text = ~pathway,
                        marker = list(
                          size = ~(-log10(pval)*8)
                        ),
                        color = ~(log10(pval)))

      container[[paste("int_plot", dataset, sep = "_")]] <- i_plot
    }

  }

  else {

    container <- huva_exp[[study]][["gsea"]]

    container_tmp <- list()

    for (i in names(container)) {

      gsea <- container[[i]]

      plot <- ggplot(gsea, aes(x = NES, y = -log10(pval), size = -log10(pval), alpha = -log10(pval)))+
        geom_point() + theme_bw() + theme(aspect.ratio = 2) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_x_continuous(limits = c(-max(abs(gsea$NES)), max(abs(gsea$NES))))

      container_tmp[[i]] <- gsea
      container_tmp[[paste("plot", i, sep = "_")]] <- plot

      if (int_plot==T) {
        i_plot <- plot_ly(type = "scatter",
                          mode = 'markers',
                          data = gsea,
                          x = ~NES, y = ~(-log10(pval)),
                          text = ~pathway,
                          marker = list(
                            size = ~(-log10(pval)*8)
                          ),
                          color = ~(log10(pval)))

        container_tmp[[paste("int_plot", i, sep = "_")]] <- i_plot
      }

    }

    container <- container_tmp

  }

  return(container)

}


#' get_anno.stats_huva
#'
#' @title get_anno.stats_huva
#' @description Differences in the annotation parameters within the two groups are statistically explored with the function get_anno.stat.huva.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @param easytoread if TRUE (default), it transforms all outputs in data frame objects.
#' @return returns a list of statistics performed on the annotation table.
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
#' anno.stat <- get_anno.stat_huva(huva_exp = binned_dataset,
#'                                 study = "FG500")
#'
#' @export
get_anno.stat_huva <- function(huva_exp, study, dataset = NULL, easytoread = T) {

  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }

  if (length(grep(study, names(huva_exp[["summary"]][["anno"]])))<1) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[["summary"]][["anno"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    container <- huva_exp[["summary"]][["anno"]][[dataset]]

    if (easytoread==TRUE) {

      for (i in names(container)) {

        if (class(container[[i]])!="htest") {

          container[[i]] <- as.data.frame(container[[i]])

        } else {

          test_table <- container[[i]]
          test_table <- as.data.frame(c(test_table$p.value, test_table$estimate))
          colnames(test_table) <- i
          rownames(test_table) <- c("p value", "mean high", "mean low") # This 4 lines are the bug fix from v 0.1.5
          test_table$stat <- rownames(test_table)
          test_table <- test_table[, c("stat", i)]
          rownames(test_table) <- NULL

          container[[i]] <- test_table


        }
      }
    }

  }
  else {
    container <- huva_exp[["summary"]][["anno"]][grep(study, names(huva_exp[["summary"]][["anno"]]))]

    if (easytoread==TRUE) {

      for (j in names(container)) {

        working_on_list <- container[[j]]

        for (i in names(working_on_list)) {

          if (class(working_on_list[[i]])!="htest") {

            working_on_list[[i]] <- as.data.frame(working_on_list[[i]])

          } else {

            test_table <- working_on_list[[i]]
            test_table <- as.data.frame(c(test_table$p.value, test_table$estimate))
            colnames(test_table) <- i
            rownames(test_table) <- c("p value", "mean high", "mean low") # This 4 lines are the bug fix from v 0.1.5
            test_table$stat <- rownames(test_table)
            test_table <- test_table[, c("stat", i)]
            rownames(test_table) <- NULL

            working_on_list[[i]] <- test_table

          }

        }

        container[[j]] <- working_on_list

      }

    }

  }

  return(container)
}

#' get_meta.stats_huva
#'
#' @title get_meta.stats_huva
#' @description Within the two huva experiment groups, statistical differences in metadata parameters can be explored with the function get_meta.stat.huva.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @return returns a list of statistics performed on the annotation table.
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
#' meta.stat <- get_meta.stat_huva(huva_exp = binned_dataset,
#'                                 study = "FG500",
#'                                 dataset = "FG500_whole_blood_cellcount")
#'
#' @export
get_meta.stat_huva <- function (huva_exp, study, dataset = NULL) {

  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }
  if (length(grep(study, names(huva_exp[["summary"]][["metadata"]]))) <
      1) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }
  if (!is.null(dataset)) {
    if (!dataset %in% names(huva_exp[["summary"]][["metadata"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }
    container <- huva_exp[["summary"]][["metadata"]][[dataset]]
  }
  else {
    container <- huva_exp[["summary"]][["metadata"]][grep(study, names(huva_exp[["summary"]][["metadata"]]))]
  }
  return(container)
}

#' get_anno.plot_huva
#'
#' @title get_anno.plot_huva
#' @description The distribution of donors in the two huva experiment groups is correlated to the available annotation parameters with the get_anno.plot_huva function.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @return The function returns a plot list based on the annotation table.
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
#' anno.plot <- get_anno.plot_huva(huva_exp = binned_dataset,
#'                                 study = "FG500")
#'
#' anno.plot$FG500_whole_blood$living_situation
#'
#' @export
get_anno.plot_huva <- function(huva_exp, study, dataset = NULL) {

  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["anno"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    data_tmp <- huva_exp[[study]][["anno"]][[dataset]]

    container <- list()

    for (w in colnames(data_tmp)[-c(1,2,3)]) {

      if (is.numeric(data_tmp[[w]])) {
        df_tmp <- data_tmp[,c("group", w)]
        colnames(df_tmp) <- c("group", "variable")
        plot <- ggplot(df_tmp, aes(x = group, y = variable, fill = group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(paste(w, study, dataset, sep = " ")) +
          theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
          scale_fill_manual(values = c("#c93918", "#0e2e99"))

      } else {
        df_tmp <- data_tmp[,c("group",w)]
        colnames(df_tmp) <- c("group", "variable")
        df_tmp <- table(df_tmp)
        df_tmp <- prop.table(df_tmp, margin = 1)*100
        df_tmp <- melt(df_tmp)
        plot <- ggplot(df_tmp, aes(x = variable, y = value, fill = group)) + geom_bar(position = "dodge", stat = "identity", color = "black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(paste(w, study, dataset, sep = " ")) +
          theme_minimal() + theme(aspect.ratio = 1) +
          scale_fill_manual(values = c("#c93918", "#0e2e99"))
      }

      container[[w]] <- plot

    }

  }
  else {
    data_tmp <- huva_exp[[study]][["anno"]]

    container <- list()

    for (n in names(data_tmp)) {

      data_tmp_sub <- data_tmp[[n]]

      for (w in colnames(data_tmp_sub)[-c(1,2,3)]) {

        if (is.numeric(data_tmp_sub[[w]])) {
          df_tmp <- data_tmp_sub[,c("group", w)]
          colnames(df_tmp) <- c("group", "variable")
          plot <- ggplot(df_tmp, aes(x = group, y = variable, fill = group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(paste(w, study, sep = " ")) +
            theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
            scale_fill_manual(values = c("#c93918", "#0e2e99"))

        } else {
          df_tmp <- data_tmp_sub[,c("group",w)]
          colnames(df_tmp) <- c("group", "variable")
          df_tmp <- table(df_tmp)
          df_tmp <- prop.table(df_tmp, margin = 1)*100
          df_tmp <- melt(df_tmp)
          plot <- ggplot(df_tmp, aes(x = variable, y = value, fill = group)) + geom_bar(position = "dodge", stat = "identity", color = "black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(paste(w, study, sep = " ")) +
            theme_minimal() + theme(aspect.ratio = 1) +
            scale_fill_manual(values = c("#c93918", "#0e2e99"))
        }

        container[[n]][[w]] <- plot

      }

    }
  }

  return(container)

}

#' get_meta.plot_huva
#'
#' @title get_meta.plot_huva
#' @description The graphical representation of the distribution of individuals across the two huva groups can be visualiuzed in correlation to metadata parameters with the get_anno.plot_huva function.
#' @param huva_exp huva experiment class object.
#' @param study character vector defining the name of the studies to be visualised.
#' @param dataset if NULL (default) all datasets will be used.
#' @return The function returns a plot list on the metadata table.
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
#' meta.plot <- get_meta.plot_huva(huva_exp = binned_dataset,
#'                                 study = "FG500")
#'
#' meta.plot$FG500_whole_blood_cellcount$Granulocytes
#'
#' @export
get_meta.plot_huva <- function(huva_exp, study, dataset = NULL) {

  if (class(huva_exp)!= "huva_experiment") {
    print("Invalid class of huva experiment used as input, it can lead to incorrect results.")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found among those provided by the huva experiment.")
    stop()
  }

  if (! "metadata" %in% names(huva_exp[[study]])) {
    error("No metadata available for this study.")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["metadata"]])) {
      print("Provided dataset not found in the specified huva experiment.")
      stop()
    }

    data_tmp <- huva_exp[[study]][["metadata"]][[dataset]]

    container <- list()

    for (w in colnames(data_tmp)[-c(1,2,3)]) {

      if (is.numeric(data_tmp[[w]])) {
        df_tmp <- data_tmp[,c("group", w)]
        colnames(df_tmp) <- c("group", "variable")
        plot <- ggplot(df_tmp, aes(x = group, y = variable, fill = group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(paste(study, dataset, sep = " ")) +
          theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
          scale_fill_manual(values = c("#c93918", "#0e2e99"))

      } else {
        df_tmp <- data_tmp[,c("group",w)]
        colnames(df_tmp) <- c("group", "variable")
        df_tmp <- table(df_tmp)
        df_tmp <- prop.table(df_tmp, margin = 1)*100
        df_tmp <- melt(df_tmp)
        plot <- ggplot(df_tmp, aes(x = variable, y = value, fill = group)) + geom_bar(position = "dodge", stat = "identity", color = "black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(paste(study, dataset, sep = " ")) +
          theme_minimal() + theme(aspect.ratio = 1) +
          scale_fill_manual(values = c("#c93918", "#0e2e99"))
      }

      container[[w]] <- plot

    }

  }
  else {
    data_tmp <- huva_exp[[study]][["metadata"]]

    container <- list()

    for (n in names(data_tmp)) {

      data_tmp_sub <- data_tmp[[n]]

      for (w in colnames(data_tmp_sub)[-c(1,2,3)]) {

        if (is.numeric(data_tmp_sub[[w]])) {
          df_tmp <- data_tmp_sub[,c("group", w)]
          colnames(df_tmp) <- c("group", "variable")
          plot <- ggplot(df_tmp, aes(x = group, y = variable, fill = group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(study) +
            theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
            scale_fill_manual(values = c("#c93918", "#0e2e99"))

        } else {
          df_tmp <- data_tmp_sub[,c("group",w)]
          colnames(df_tmp) <- c("group", "variable")
          df_tmp <- table(df_tmp)
          df_tmp <- prop.table(df_tmp, margin = 1)*100
          df_tmp <- melt(df_tmp)
          plot <- ggplot(df_tmp, aes(x = variable, y = value, fill = group)) + geom_bar(position = "dodge", stat = "identity", color = "black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(study) +
            theme_minimal() + theme(aspect.ratio = 1) +
            scale_fill_manual(values = c("#c93918", "#0e2e99"))
        }

        container[[n]][[w]] <- plot

      }

    }
  }

  return(container)

}
