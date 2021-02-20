#' Visualization of the expression histogram for a selected gene of interest
#'
#' @title get_expr.plot_exam
#' @description The function plots expression of the gene of interest (GOI) across all datasets to inspect possible
#'              distribution differences in a huva gene examination experiment.
#' @param huva_expression huva_expression class object
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @param bins graphical option, see ggplot2 for details (default = 30).
#' @param alpha graphical option, see ggplot2 for details (default = 1).
#' @return ggplot class object, histogram reporting the expression of a selected GOI across selected datasets.
#' @seealso gene_exam, run_huva_experiment
#' @import ggplot2
#' @import reshape2
#' @import ggsci
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' expr_exam.plot <- get_expr.plot_exam(huva_expression = gene_overview)
#'
#' @export
get_expr.plot_exam <- function(huva_expression, study = "ALL", dataset = "ALL", bins = 30, alpha = 1) {

  if (class(huva_expression)!= "huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results.")
  }

  name <- huva_expression$gene_name

  table <- huva_expression[["expression"]][["expression_summary"]]
  table <- suppressMessages(melt(table, variable.name = "dataset", value.name = "expression"))
  table <- table[!is.na(table$expression),]
  table$study <- unlist(lapply(strsplit(as.character(paste(table$dataset)), split = "_"), `[[`, 1))

  if (paste(study, collapse = "")=="ALL") {

    print("plotting expression for all available datasets")

    ggplot(table, aes(x = expression, color = study))+
      geom_histogram(aes(y = ..density..),bins = bins, alpha = alpha, fill = "white") +
      geom_density(alpha = .2, color = "black", aes(fill = study))+
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

      ggplot(table, aes(x = expression, color = study))+
        geom_histogram(aes(y = ..density..), bins = bins, alpha = alpha, fill = "white") +
        geom_density(alpha = .2, color = "black", aes(fill = study))+
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

      ggplot(table, aes(x = expression, color = study))+
        geom_histogram(aes(y = ..density..), bins = bins, alpha = alpha, fill = "white") +
        geom_density(alpha = .2, color = "black", aes(fill = study))+
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
#' @title get_expr_exam
#' @description With this function, the gene of interest (GOI) expression table in huva gene examination selected
#'              datasets is exported as a data.frame reporting sample names as rownames and expression values in
#'              a single column.
#' @param huva_expression huva_expression class object.
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @return data.frame or list of the expression of a selected gene of interest across the available datasets.
#' @seealso gene_exam, run_huva_experiment
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' expr_exam <- get_expr_exam(huva_expression = gene_overview,
#'                            study = "ImmVar",
#'                            dataset = "ImmVar_CD4T")
#'
#' @export
get_expr_exam <- function(huva_expression, study, dataset) {

  if (class(huva_expression)!= "huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results.")
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
#' @title get_anno_exam
#' @description This function exports the comprehensive table reporting all the available annotation parameters
#'              for the selected dataset in a huva gene examination experiment.
#' @param huva_expression huva_expression class object.
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets.
#' @seealso gene_exam, run_huva_experiment
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' anno_exam <- get_anno_exam(huva_expression = gene_overview,
#'                            study = "CEDAR",
#'                            dataset = "CEDAR_CD4T")
#'
#' @export
get_anno_exam <- function(huva_expression, study = "ALL", dataset = "ALL") {

  if (class(huva_expression)!= "huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results.")
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
#' @title get_anno.stat_exam
#' @description The get_anno.stat_exam function performs a statistical correlation between the GOI expression and
#'              each of the parameters provided in the annotation table of a huva gene examination experiment .
#' @param huva_expression huva_expression class object.
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @param param colName of the sample table to be used for the statistical analysis.
#' @return list of statistics on the annotation data of the huva examination experiment.
#' @seealso gene_exam, run_huva_experiment
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' anno.stat_exam <- get_anno.stat_exam(huva_expression = gene_overview,
#'                                      study = "FG500", dataset = "ALL")
#'
#' @export
get_anno.stat_exam <- function(huva_expression, study = "ALL", dataset = "ALL", param = "ALL"){

  if (class(huva_expression)!= "huva_gene_examination") {
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

        if (length(levels(tmp_df$param))>2 & length(levels(tmp_df$param))!= length(tmp_df$param)){

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
#' @title get_anno.plot_exam
#' @description Correlations in a huva gene examination experiment can be graphically visualized with the
#'              get_anno.plot_exam function, whose format produces an output containing a list of graphical objects.
#' @param huva_expression huva_expression class object.
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @return The function returns a plot list based on the annotation table.
#' @seealso gene_exam, run_huva_experiment
#' @import ggplot2
#' @import ggpubr
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' anno.plot_exam <- get_anno.plot_exam(huva_expression = gene_overview,
#'                                      study = "FG500",
#'                                      dataset = "FG500_whole_blood")
#'
#' anno.plot_exam$weight
#'
#' @export
get_anno.plot_exam <- function(huva_expression, study = "ALL", dataset = "ALL"){

  if (class(huva_expression)!= "huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results.")
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

          plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
            geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

        }

        else {

          tmp <- data_anno

          colnames(tmp)[colnames(tmp)==j] <- "param"
          plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
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

            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_anno

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
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

            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_anno

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
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
#' @title get_meta_exam
#' @description The function provides metadata table parameters for a selected gene of interest used in a huva
#'              gene examination experiment.
#' @param huva_expression huva_expression class object.
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @return data.frame or list of the annotation data of a selected gene of interest across the available datasets.
#' @seealso gene_exam, run_huva_experiment
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db,
#'                            gene = "MYD88")
#'
#' meta.table_exam <- get_meta_exam(huva_expression = gene_overview,
#'                                  study = "FG500")
#'
#' @export
get_meta_exam <- function(huva_expression, study = "ALL", dataset = "ALL") {

  if (class(huva_expression)!= "huva_gene_examination") {
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
#' @title get_meta.stat_exam
#' @description Huva examination experiment statistical analysis of the influence of provided metadata on the expression of a gene of
#'              interest (GOI) is performed by the function get_meta.stat_exam which reports it on a data.frame or list of data.frames.
#' @param huva_expression huva_expression class object
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @param param colName of the sample table to be used for the statistical analysis.
#' @return list of statistics on the metadata of the huva examination experiment.
#' @seealso gene_exam, run_huva_experiment
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' meta.stat_exam <- get_meta.stat_exam(huva_expression = gene_overview,
#'                                      study = "FG500",
#'                                      dataset = "FG500_whole_blood_cellcount")
#'
#' @export
get_meta.stat_exam <- function(huva_expression, study = "ALL", dataset = "ALL", param = "ALL"){

  if (class(huva_expression)!= "huva_gene_examination") {
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
#' @title get_meta.plot_exam
#' @description The get_meta.plot_exam function will quickly provide a list of metadata-based graphical objects in
#'              a huva examination experiment.
#' @param huva_expression huva_expression class object.
#' @param study character vector defining the studies to be used in the analysis (default = "ALL").
#' @param dataset character vector defining the datasets to be used in the analysis (default = "ALL").
#' @return The function returns a plot list based on the annotation table.
#' @seealso gene_exam, run_huva_experiment
#' @import ggplot2
#' @import ggpubr
#' @examples
#' library(huva)
#' library(huva.db)
#'
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' meta.plot_exam <- get_meta.plot_exam(huva_expression = gene_overview,
#'                                      study = "FG500",
#'                                      dataset = "FG500_whole_blood_cellcount")
#'
#' meta.plot_exam$`Monocytes (CD14+)`
#'
#' @export
get_meta.plot_exam <- function(huva_expression, study = "ALL", dataset = "ALL"){

  if (class(huva_expression)!= "huva_gene_examination") {
    error("Error: Use huva_gene_examination class object for reliable results.")
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
          plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
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

            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_metadata

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
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

            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
              geom_point() + theme_minimal() + theme(aspect.ratio = 1) + ggtitle(paste(huva_expression$gene_name)) + stat_cor(method = "pearson")

          }

          else {

            tmp <- data_metadata

            colnames(tmp)[colnames(tmp)==j] <- "param"
            plot <- ggplot(tmp, aes(x = param, y = expression))+ xlab(j) +
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
