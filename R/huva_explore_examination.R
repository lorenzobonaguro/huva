#' Visualization of the expression histogram for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @param bins see ggplot2 documentation (default = 30)
#' @param alpha see ggplot2 documentation (default = 1)
#' @return ggplot class object, histogram for the expression of a selected gene of interest across selected datasets
#' @examples
#'
#' 2+2
#'
#' @import ggplot2
#' @import reshape2
#' @import ggsci
#' @export
get_expr.plot_exam <- function(huva_expression, study = "ALL", dataset = "ALL", bins=30, alpha=1) {

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  name <- huva_expression$gene_name

  table <- huva_expression[["expression"]][["expression_summary"]]
  table <- suppressMessages(melt(table, variable.name = "dataset", value.name = "expression"))
  table <- table[!is.na(table$expression),]
  table$study <- unlist(lapply(strsplit(as.character(paste(table$dataset)), split = "_"), `[[`, 1))

  if (paste(study, collapse = "")=="ALL") {

    print("plotting expression for all available datasets")

    ggplot(table, aes(x=expression, color=study))+
      geom_histogram(aes(y=..density..),bins = bins, alpha=alpha, fill="white") +
      geom_density(alpha=.2, color="black", aes(fill=study))+
      ggtitle(paste(name, "expression"))+
      theme_minimal() +
      scale_fill_aaas()+
      theme(aspect.ratio = 1)+
      scale_color_aaas()+
      facet_wrap(~dataset, scales = "free")
  }

  else {

    if (paste(dataset, collapse = "")=="ALL") {

      print(paste("plotting expression for", paste(study, collapse = ", ")))

      table <- table[table$study %in% study,]

      ggplot(table, aes(x=expression, color=study))+
        geom_histogram(aes(y=..density..), bins = bins, alpha=alpha, fill="white") +
        geom_density(alpha=.2, color="black", aes(fill=study))+
        ggtitle(paste(name, "expression"))+
        theme_minimal() +
        scale_fill_aaas()+
        theme(aspect.ratio = 1)+
        scale_color_aaas()+
        facet_wrap(~dataset, scales = "free")

    }

    else {

      print(paste("plotting expression for", paste(dataset, collapse = ", ")))

      table <- table[table$dataset %in% dataset,]

      ggplot(table, aes(x=expression, color=study))+
        geom_histogram(aes(y=..density..), bins = bins, alpha=alpha, fill="white") +
        geom_density(alpha=.2, color="black", aes(fill=study))+
        ggtitle(paste(name, "expression"))+
        theme_minimal() +
        scale_fill_aaas()+
        theme(aspect.ratio = 1)+
        scale_color_aaas()+
        facet_wrap(~dataset, scales = "free")

    }

  }

}

#' Expression table for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @return data.frame or list of the expression of a selected gene of interest across the available datasets
#' @examples
#'
#' 2+2
#'
#' @export
get_expr_exam <- function(huva_expression, study, dataset) {

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  if (paste(study, collapse = "")=="ALL") {

    table <- huva_expression[["expression"]][["expression_summary"]]

  } else{

    table <- huva_expression[["expression"]][grep(names(huva_expression[["expression"]]), pattern = study)]

    if (paste(dataset, collapse = "")!="ALL") {

      table <- table[dataset]

    }

  }

  if (length(table)==1) {
    table <- as.data.frame(table)
  }

  return(table)

}

#' Annotation table for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets
#' @examples
#'
#' 2+2
#'
#' @export
get_anno_exam <- function(huva_expression, study="ALL", dataset="ALL") {

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  if (paste(study, collapse = "")=="ALL") {
    table <- huva_expression[["anno"]]
  } else {

    table <- huva_expression[["anno"]][grep(names(huva_expression[["anno"]]), pattern = study)]

    if (paste(dataset, collapse = "")!="ALL") {

      table <- table[dataset]
    }

  }

  if (length(table)==1) {
    table <- as.data.frame(table)
    colnames(table) <- lapply(strsplit(colnames(table), split="\\."), `[`, 2)
    colnames(table)[1] <- "Row.names"
  }

  return(table)

}

#' Statistical analysis of annotation table for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @param param colName of the sample table to be used for the statistical analyisis
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets
#' @examples
#'
#' 2+2
#'
#' @export
get_anno.stat_exam <- function(huva_expression, study="ALL", dataset="ALL", param="ALL"){

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  if (paste(study, collapse = "")=="ALL") {

    container <- list()

  }

  else {
    if (paste(dataset, collapse = "")=="ALL") {

      huva_expression[["anno"]] <- huva_expression[["anno"]][grep(names(huva_expression[["anno"]]), pattern = study)]

      container <- list()

    } else {

      huva_expression[["anno"]] <- huva_expression[["anno"]][dataset]

      container <- list()

    }
  }

  for (i in names(huva_expression[["anno"]])) {

    data_anno <- huva_expression[["anno"]][[i]]

    df_container <- matrix(nrow = 3)
    df_container <- as.data.frame(df_container)

    for (j in colnames(data_anno)[-c(1,2)]) {

      if (is.numeric(data_anno[[j]])) {

        p.val <- cor.test(data_anno[["expression"]], data_anno[[j]])$p.value
        type <- "pearson"
        vector <- c(j, p.val, type)

      } else {

        tmp_df <- data_anno
        colnames(tmp_df)[colnames(tmp_df)==j] <- "param"
        tmp_df$param <- as.factor(tmp_df$param)

        if(length(levels(tmp_df$param))==1) {

          vector <- c(j, "ND", "only one variable")

        }

        if (length(levels(tmp_df$param))==2) {

          p.val <- t.test(expression~param, data = tmp_df)$p.value
          type <- "t.test"
          vector <- c(j, p.val, type)

        }

        if(length(levels(tmp_df$param))==length(tmp_df$param)){

          paste("get")

          vector <- c(j, "ND", "all unique variables")

        }

        if (length(levels(tmp_df$param))>2 & length(levels(tmp_df$param))!=length(tmp_df$param)){

          p.val <- summary(aov(expression~param, data = tmp_df))[[1]][["Pr(>F)"]][1]
          type <- "one-way-anova"
          vector <- c(j, p.val, type)

        }

      }

      df_container[[j]] <- vector

    }

    df_container$V1 <- NULL
    df_container <- t(df_container)
    rownames(df_container) <- NULL
    colnames(df_container) <- c("Param", "p.val", "Test")

    container[[i]] <- df_container

  }

  return(container)

}

#' Plot of annotation table parameters for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets
#' @examples
#'
#' 2+2
#'
#' @import ggplot2
#' @import ggpubr
#' @export
get_anno.plot_exam <- function(huva_expression, study="ALL", dataset="ALL"){

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  if (paste(study, collapse = "")=="ALL") {

    container <- list()

    for (i in names(huva_expression[["anno"]])) {

      data_anno <- huva_expression[["anno"]][[i]]

      container_plot <- list()

      for (j in colnames(data_anno)[-c(1,2)]) {

        if(is.numeric(data_anno[[j]])) {

          tmp <- data_anno

          colnames(tmp)[colnames(tmp)==j] <- "param"

          plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
            geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

        }

        else {

          tmp <- data_anno

          colnames(tmp)[colnames(tmp)==j] <- "param"
          plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
            geom_boxplot()+ theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name))

        }

        container_plot[[j]] <- plot

      }

      container[[i]] <- container_plot

    }

  }

  else {

    if (paste(dataset, collapse = "")=="ALL") {

      huva_expression[["anno"]] <- huva_expression[["anno"]][grep(names(huva_expression[["anno"]]), pattern = study)]

      container <- list()

      for (i in names(huva_expression[["anno"]])) {

        data_anno <- huva_expression[["anno"]][[i]]

        container_plot <- list()

        for (j in colnames(data_anno)[-c(1,2)]) {

          if(is.numeric(data_anno[[j]])) {

            tmp <- data_anno

            colnames(tmp)[colnames(tmp)==j] <- "param"

            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_anno

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_boxplot()+ theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name))

          }

          container_plot[[j]] <- plot

        }

        container[[i]] <- container_plot

      }

    }

    else {

      huva_expression[["anno"]] <- huva_expression[["anno"]][dataset]

      container <- list()

      for (i in names(huva_expression[["anno"]])) {

        data_anno <- huva_expression[["anno"]][[i]]

        container_plot <- list()

        for (j in colnames(data_anno)[-c(1,2)]) {

          if(is.numeric(data_anno[[j]])) {

            tmp <- data_anno

            colnames(tmp)[colnames(tmp)==j] <- "param"

            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_anno

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_boxplot()+ theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name))

          }

          container_plot[[j]] <- plot

        }

        container <- container_plot

      }

    }

  }

  # here add the other if statement

  return(container)

}

#' Metadata table parameters for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets
#' @examples
#'
#' 2+2
#'
#' @export
get_meta_exam <- function(huva_expression, study="ALL", dataset="ALL") {

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  if (paste(study, collapse = "")=="ALL") {
    table <- huva_expression[["metadata"]]
  } else {

    table <- huva_expression[["metadata"]][grep(names(huva_expression[["metadata"]]), pattern = study)]

    if (paste(dataset, collapse = "")!="ALL") {

      table <- table[dataset]
    }

  }

  if (length(table)==1) {
    table <- as.data.frame(table)
    colnames(table) <- lapply(strsplit(colnames(table), split="\\."), `[`, 2)
    colnames(table)[1] <- "Row.names"
  }

  return(table)

}

#' Metadata table parameters for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @param param colName of the sample table to be used for the statistical analyisis
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets
#' @examples
#'
#' 2+2
#'
#' @export
get_meta.stat_exam <- function(huva_expression, study="ALL", dataset="ALL", param="ALL"){

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  if (paste(study, collapse = "")=="ALL") {

    container <- list()

  }

  else {
    if (paste(dataset, collapse = "")=="ALL") {

      huva_expression[["metadata"]] <- huva_expression[["metadata"]][grep(names(huva_expression[["metadata"]]), pattern = study)]

      container <- list()

    } else {

      huva_expression[["metadata"]] <- huva_expression[["metadata"]][dataset]

      container <- list()

    }
  }

  for (i in names(huva_expression[["metadata"]])) {

    data_metadata <- huva_expression[["metadata"]][[i]]

    df_container <- matrix(nrow = 3)
    df_container <- as.data.frame(df_container)

    for (j in colnames(data_metadata)[-c(1,2)]) {

      data_metadata_tmp <- data_metadata[!is.na(data_metadata[[j]]),]

      if (is.numeric(data_metadata_tmp[[j]])) {

        p.val <- cor.test(data_metadata_tmp[["expression"]], data_metadata_tmp[[j]])$p.value
        type <- "pearson"
        vector <- c(j, p.val, type)

      } else {

        tmp_df <- data_metadata_tmp
        colnames(tmp_df)[colnames(tmp_df)==j] <- "param"
        tmp_df$param <- as.factor(tmp_df$param)

        if(length(levels(tmp_df$param))==1) {

          vector <- c(j, "ND", "only one variable")

        }

        if (length(levels(tmp_df$param))==2) {

          p.val <- t.test(expression~param, data = tmp_df)$p.value
          type <- "t.test"
          vector <- c(j, p.val, type)

        }

        if(length(levels(tmp_df$param))==length(tmp_df$param)){

          paste("all unique")

          vector <- c(j, "ND", "all unique variables")

        }

        if (length(levels(tmp_df$param))>2 & length(levels(tmp_df$param))!=length(tmp_df$param)){

          p.val <- summary(aov(expression~param, data = tmp_df))[[1]][["Pr(>F)"]][1]
          type <- "one-way-anova"
          vector <- c(j, p.val, type)

        }

      }

      df_container[[j]] <- vector

    }

    df_container$V1 <- NULL
    df_container <- t(df_container)
    rownames(df_container) <- NULL
    colnames(df_container) <- c("Param", "p.val", "Test")

    container[[i]] <- df_container

  }

  return(container)

}

#' Metadata table parameters for a selected gene of interest
#'
#' @param huva_expression huva_expression class object
#' @param study charachter or charachter vector defining the studies to be used in the analysis (default = "ALL")
#' @param dataset charachter or charachter vector defining the datasets to be used in the analysis (default = "ALL")
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets
#' @examples
#'
#' 2+2
#'
#' @import ggplot2
#' @import ggpubr
#' @export
get_meta.plot_exam <- function(huva_expression, study="ALL", dataset="ALL"){

  if (class(huva_expression)!="huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results")
  }

  if (paste(study, collapse = "")=="ALL") {

    container <- list()

    for (i in names(huva_expression[["metadata"]])) {

      data_metadata <- huva_expression[["metadata"]][[i]]

      container_plot <- list()

      for (j in colnames(data_metadata)[-c(1,2)]) {

        if(is.numeric(data_metadata[[j]])) {

          tmp <- data_metadata

          colnames(tmp)[colnames(tmp)==j] <- "param"

          plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
            geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

        }

        else {

          tmp <- data_metadata

          colnames(tmp)[colnames(tmp)==j] <- "param"
          plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
            geom_boxplot()+ theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name))

        }

        container_plot[[j]] <- plot

      }

      container[[i]] <- container_plot

    }

  }

  else {

    if (paste(dataset, collapse = "")=="ALL") {

      huva_expression[["metadata"]] <- huva_expression[["metadata"]][grep(names(huva_expression[["metadata"]]), pattern = study)]

      container <- list()

      for (i in names(huva_expression[["metadata"]])) {

        data_metadata <- huva_expression[["metadata"]][[i]]

        container_plot <- list()

        for (j in colnames(data_metadata)[-c(1,2)]) {

          if(is.numeric(data_metadata[[j]])) {

            tmp <- data_metadata

            colnames(tmp)[colnames(tmp)==j] <- "param"

            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_metadata

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_boxplot()+ theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name))

          }

          container_plot[[j]] <- plot

        }

        container[[i]] <- container_plot

      }

    }

    else {

      huva_expression[["metadata"]] <- huva_expression[["metadata"]][dataset]

      container <- list()

      for (i in names(huva_expression[["metadata"]])) {

        data_metadata <- huva_expression[["metadata"]][[i]]

        container_plot <- list()

        for (j in colnames(data_metadata)[-c(1,2)]) {

          if(is.numeric(data_metadata[[j]])) {

            tmp <- data_metadata

            colnames(tmp)[colnames(tmp)==j] <- "param"

            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_metadata

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x=param, y=expression))+ xlab(j) +
              geom_boxplot()+ theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name))

          }

          container_plot[[j]] <- plot

        }

        container <- container_plot

      }

    }

  }

  # here add the other if statement

  return(container)

}
