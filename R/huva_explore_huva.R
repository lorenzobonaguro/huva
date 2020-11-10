#' get_expr_huva
#'
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @return expression table of the binned samples
#' @examples
#'
#' 2+2
#'
#' @export
get_expr_huva <- function(huva_exp, study, dataset=NULL) {

  if (class(huva_exp)!="huva_experiment") {
    error("Worng class of huva experiment used. Wrong results may occur")
  }

  if (! study %in% names(huva_exp)) {
    error("Study not found in the huva experiment provided")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["data"]])) {
      error("Dataset not found in the huva experiment provided")
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
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @return annotation table of the binned samples
#' @examples
#'
#' 2+2
#'
#' @export
get_anno_huva <- function(huva_exp, study, dataset=NULL) {

  if (class(huva_exp)!="huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found in the huva experiment provided")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["anno"]])) {
      print("Dataset not found in the huva experiment provided")
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
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @return annotation table of the binned samples
#' @examples
#'
#' 2+2
#'
#' @export
get_meta_huva <- function(huva_exp, study, dataset=NULL) {

  if (class(huva_exp)!="huva_experiment") {
    error("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    error("Study not found in the huva experiment provided")
    stop()
  }

  if (! "metadata" %in% names(huva_exp[[study]])) {
    error("No metadata avaiilable for this study")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["metadata"]])) {
      print("Dataset not found in the huva experiment provided")
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
#' @param huva_exp huva experiment class object
#' @param study character vector defining the names of the study to be visualized
#' @param dataset if NULL (default) all datasets will be used
#' @param cluster_col TRUE defaul
#' @param pval p-value cut-off used for the definition of differential expression
#' @param logFC log2 fold change cut-off for the definition of differential expression
#' @param PC character vetor defining the PC to be plotted (c("PC1", "PC2"), default)
#' @return table of deferentially expressed genes in the comparison between the binned samples
#' @import pheatmap
#' @import grDevices
#' @examples
#'
#' 2+2
#'
#' @export
get_DE_huva <- function(huva_exp, study, dataset=NULL, pval=0.05, logFC=0,cluster_col=T, PC = c("PC1", "PC2")) {

  if (class(huva_exp)!="huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found in the huva experiment provided")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["DE_genes"]])) {
      print("Dataset not found in the huva experiment provided")
      stop()
    }

    container <- list()
    data <- huva_exp[[study]][["DE_genes"]][[dataset]]
    data <- data[data$adj.P.Val < pval, ]
    data <- data[abs(data$logFC) > logFC, ]
    data$direction <- ifelse(data$logFC < 0, "down", "up")
    container[[dataset]] <- data

    container[[paste("plot", dataset, sep = "_")]] <- ggplot(as.data.frame(table(data$direction)), aes(x=Var1, y=Freq, fill=Var1))+
      geom_bar(stat = "identity", colour="black", size=0.5) + theme_minimal() +
      geom_text(stat = "identity", aes(label=Freq), vjust= -0.2) +
      theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_text(size = 12), axis.title = element_text(size=14)) +
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
                              annotation_colors = list(group=c(high="#c93918", low="#0e2e99"),
                                                       expression=heat.colors(20, alpha = 1)), silent = T)

    # PCA
    pca_data <- huva_exp[[study]][["data"]][[dataset]]
    pca_meta <- hm_meta

    pc <- prcomp(t(pca_data))$x
    pc <- merge(pc, pca_meta, by="row.names")

    container[[paste("PCA", dataset, sep = "_")]] <- ggplot(pc, aes(x=.data[[PC[1]]], y=.data[[PC[2]]], colour=group))+
                                                        geom_point(size=5) +
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
      container_tmp[[paste("plot", i, sep = "_")]] <- ggplot(as.data.frame(table(data$direction)), aes(x=Var1, y=Freq, fill=Var1))+
        geom_bar(stat = "identity", colour="black", size=0.5) + theme_minimal() +
        geom_text(stat = "identity", aes(label=Freq), vjust= -0.2) +
        theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_text(size = 12), axis.title = element_text(size=14)) +
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
                                                             annotation_colors = list(group=c(high="#c93918", low="#0e2e99"),
                                                                                      expression=heat.colors(20, alpha = 1)), silent = T)

      pca_data <- container[["data"]][[i]]
      pca_meta <- hm_meta

      pc <- prcomp(t(pca_data))$x
      pc <- merge(pc, pca_meta, by="row.names")

      container_tmp[[paste("PCA", i, sep = "_")]] <- ggplot(pc, aes(x=.data[[PC[1]]], y=.data[[PC[2]]], colour=group))+
        geom_point(size=5) +
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
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @param n_top_genes number of top down or up regulated genes to be displayed
#' @return ranking of the genes
#' @examples
#'
#' 2+2
#'
#' @export
get_rank_huva <- function(huva_exp, study, dataset=NULL, n_top_genes=20) {

  if (class(huva_exp)!="huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found in the huva experiment provided")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["Rank_genelist"]])) {
      print("Dataset not found in the huva experiment provided")
      stop()
    }

    rank <- huva_exp[[study]][["Rank_genelist"]][[dataset]]
    rank_top <- rbind(data.frame(delta_exp = rank[1:n_top_genes], gene_name = names(rank[1:n_top_genes]), group = "UP"),
                      data.frame(delta_exp = rank[(length(rank)-n_top_genes+1):length(rank)], gene_name = names(rank[(length(rank)-n_top_genes+1):length(rank)]), group = "DOWN"))

    rank_top$gene_name <- factor(rank_top$gene_name, levels = rank_top$gene_name)

    plot <- ggplot(rank_top, aes(x=gene_name, y=delta_exp, fill=group)) +
      geom_point(size=6, shape=21) + facet_wrap(~group, scales = "free") + coord_flip() +
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

      plot <- ggplot(rank_top, aes(x=gene_name, y=delta_exp, fill=group)) +
        geom_point(size=6, shape=21) + facet_wrap(~group, scales = "free") + coord_flip() +
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
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @param int_plot default TRUE, decide if the interactive GSEA plot will be calculated
#' @return return GSEA table and plot
#' @examples
#'
#' 2+2
#'
#' @importFrom  plotly plot_ly
#' @export
get_gsea_huva <- function(huva_exp, study, dataset=NULL, int_plot=T){

  if (class(huva_exp)!="huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found in the huva experiment provided")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["gsea"]])) {
      print("Dataset not found in the huva experiment provided")
      stop()
    }

    gsea <- huva_exp[[study]][["gsea"]][[dataset]]

    plot <- ggplot(gsea, aes(x=NES, y=-log10(pval), size=-log10(pval), alpha=-log10(pval)))+
      geom_point() + theme_bw() + theme(aspect.ratio = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_x_continuous(limits = c(-max(abs(gsea$NES)), max(abs(gsea$NES))))

    container <- list()

    container[[dataset]] <- gsea
    container[[paste("plot", dataset, sep = "_")]] <- plot

    if (int_plot==T) {
      i_plot <- plot_ly(type="scatter",
                        mode = 'markers',
                        data = gsea,
                        x=~NES, y=~(-log10(pval)),
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

      plot <- ggplot(gsea, aes(x=NES, y=-log10(pval), size=-log10(pval), alpha=-log10(pval)))+
        geom_point() + theme_bw() + theme(aspect.ratio = 2) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_x_continuous(limits = c(-max(abs(gsea$NES)), max(abs(gsea$NES))))

      container_tmp[[i]] <- gsea
      container_tmp[[paste("plot", i, sep = "_")]] <- plot

      if (int_plot==T) {
        i_plot <- plot_ly(type="scatter",
                          mode = 'markers',
                          data = gsea,
                          x=~NES, y=~(-log10(pval)),
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
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @param easytoread if TRUE (default) transforms all output in data.frame ojects
#' @return return statistics list on annotation table
#' @examples
#'
#' 2+2
#'
#' @export
get_anno.stat_huva <- function(huva_exp, study, dataset=NULL, easytoread=T) {

  if (class(huva_exp)!="huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (length(grep(study, names(huva_exp[["summary"]][["anno"]])))<1) {
    print("Study not found in the huva experiment provided")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[["summary"]][["anno"]])) {
      print("Dataset not found in the huva experiment provided")
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
          rownames(test_table)[1] <- "p value"

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
            rownames(test_table)[1] <- "p value"

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
#' @param huva_exp huva experiment class object
#' @param study Character vector defining the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @return return statistics list on metadata table
#' @examples
#'
#' 2+2
#'
#' @export
get_meta.stat_huva <- function (huva_exp, study, dataset = NULL) {

  if (class(huva_exp) != "huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }
  if (length(grep(study, names(huva_exp[["summary"]][["metadata"]]))) <
      1) {
    print("Study not found in the huva experiment provided")
    stop()
  }
  if (!is.null(dataset)) {
    if (!dataset %in% names(huva_exp[["summary"]][["metadata"]])) {
      print("Dataset not found in the huva experiment provided")
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
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @return return plot list on annotation table
#' @examples
#'
#' 2+2
#'
#' @export
get_anno.plot_huva <- function(huva_exp, study, dataset=NULL) {

  if (class(huva_exp)!="huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found in the huva experiment provided")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["anno"]])) {
      print("Dataset not found in the huva experiment provided")
      stop()
    }

    data_tmp <- huva_exp[[study]][["anno"]][[dataset]]

    container <- list()

    for (w in colnames(data_tmp)[-c(1,2,3)]) {

      if (is.numeric(data_tmp[[w]])) {
        df_tmp <- data_tmp[,c("group", w)]
        colnames(df_tmp) <- c("group", "variable")
        plot <- ggplot(df_tmp, aes(x=group, y=variable, fill=group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(paste(w, study, dataset, sep = " ")) +
          theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
          scale_fill_manual(values = c("#c93918", "#0e2e99"))

      } else {
        df_tmp <- data_tmp[,c("group",w)]
        colnames(df_tmp) <- c("group", "variable")
        df_tmp <- table(df_tmp)
        df_tmp <- prop.table(df_tmp, margin = 1)*100
        df_tmp <- melt(df_tmp)
        plot <- ggplot(df_tmp, aes(x=variable, y=value, fill=group)) + geom_bar(position="dodge", stat = "identity", color="black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(paste(w, study, dataset, sep = " ")) +
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
          plot <- ggplot(df_tmp, aes(x=group, y=variable, fill=group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(paste(w, study, sep = " ")) +
            theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
            scale_fill_manual(values = c("#c93918", "#0e2e99"))

        } else {
          df_tmp <- data_tmp_sub[,c("group",w)]
          colnames(df_tmp) <- c("group", "variable")
          df_tmp <- table(df_tmp)
          df_tmp <- prop.table(df_tmp, margin = 1)*100
          df_tmp <- melt(df_tmp)
          plot <- ggplot(df_tmp, aes(x=variable, y=value, fill=group)) + geom_bar(position="dodge", stat = "identity", color="black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(paste(w, study, sep = " ")) +
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
#' @param huva_exp huva experiment class object
#' @param study Chaarachter vector defininf the names of the study to be visualised
#' @param dataset if NULL (default) all datasets will be used
#' @return return plot list on metadata table
#' @examples
#'
#' 2+2
#'
#' @export
get_meta.plot_huva <- function(huva_exp, study, dataset=NULL) {

  if (class(huva_exp)!="huva_experiment") {
    print("Worng class of huva experiment used. Wrong results may occur")
    stop()
  }

  if (! study %in% names(huva_exp)) {
    print("Study not found in the huva experiment provided")
    stop()
  }

  if (! "metadata" %in% names(huva_exp[[study]])) {
    error("No metadata avaiilable for this study")
    stop()
  }

  if (!is.null(dataset)) {

    if (! dataset %in% names(huva_exp[[study]][["metadata"]])) {
      print("Dataset not found in the huva experiment provided")
      stop()
    }

    data_tmp <- huva_exp[[study]][["metadata"]][[dataset]]

    container <- list()

    for (w in colnames(data_tmp)[-c(1,2,3)]) {

      if (is.numeric(data_tmp[[w]])) {
        df_tmp <- data_tmp[,c("group", w)]
        colnames(df_tmp) <- c("group", "variable")
        plot <- ggplot(df_tmp, aes(x=group, y=variable, fill=group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(paste(study, dataset, sep = " ")) +
          theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
          scale_fill_manual(values = c("#c93918", "#0e2e99"))

      } else {
        df_tmp <- data_tmp[,c("group",w)]
        colnames(df_tmp) <- c("group", "variable")
        df_tmp <- table(df_tmp)
        df_tmp <- prop.table(df_tmp, margin = 1)*100
        df_tmp <- melt(df_tmp)
        plot <- ggplot(df_tmp, aes(x=variable, y=value, fill=group)) + geom_bar(position="dodge", stat = "identity", color="black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(paste(study, dataset, sep = " ")) +
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
          plot <- ggplot(df_tmp, aes(x=group, y=variable, fill=group))+ geom_boxplot() + xlab("") + ylab(w) + ggtitle(study) +
            theme_minimal() + theme(aspect.ratio = 1) + stat_compare_means(method = "t.test") +
            scale_fill_manual(values = c("#c93918", "#0e2e99"))

        } else {
          df_tmp <- data_tmp_sub[,c("group",w)]
          colnames(df_tmp) <- c("group", "variable")
          df_tmp <- table(df_tmp)
          df_tmp <- prop.table(df_tmp, margin = 1)*100
          df_tmp <- melt(df_tmp)
          plot <- ggplot(df_tmp, aes(x=variable, y=value, fill=group)) + geom_bar(position="dodge", stat = "identity", color="black") + xlab("") + ylab(paste(w,"%", sep = " ")) + ggtitle(study) +
            theme_minimal() + theme(aspect.ratio = 1) +
            scale_fill_manual(values = c("#c93918", "#0e2e99"))
        }

        container[[n]][[w]] <- plot

      }

    }
  }

  return(container)

}
