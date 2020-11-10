#' Calculate huva experiment on the expression variance of a selected gene of interest
#'
#' @param data huva_dataset class object.
#' @param gene gene name used for the analysis
#' @param quantiles definition of the quantile of segregation of the samples, quantiles are always simmetrical
#'     between high and low groups (e.g. \code{quantile 0.1} will use the 10th and 90th percentiles)
#' @param gs_list class list object defining gene sets to be included in the analysis (to generate this file see
#'     the documentation of fgse)
#' @param summ bolean. define if the summary of the huva experiment will be calculated
#' @param datasets_list characher vector used to filter the dataset in the data objet for the analysis
#' @param adjust.method p value adjustment method used to corret the DE genes analysis
#' @return huva_experiment
#' @examples
#'
#' 2+2
#'
#' @import limma
#' @import fgsea
#' @import Rmisc
#' @import ggpubr
#' @import reshape2
#' @import huva.db
#' @export
run_huva_experiment <- function(data=datasets, gene, quantiles, gs_list,summ=T, datasets_list=NULL, adjust.method="none") {

  if (class(data)!="huva_dataset") {
    error("Use huva_dataset class object to run the huva_experiment function")
  }

  container <- list()

  print(paste("Binning on ", gene, " expression", sep = ""))

  if (is.null(datasets_list)==F) {
    # This is to allow the selection of the datasets to use in the analysis
    data <- data[datasets_list]
  }

  for (i in names(data)) {

    for (j in names(data[[i]][["data"]])) {

      if (gene %in% rownames(data[[i]][["data"]][[j]])) {

        expr <- as.data.frame(data[[i]][["data"]][[j]][gene,])
        colnames(expr) <- c("expression")

        # Calculation of the percentiles

        expr$group <- ifelse(expr$expression<= quantile(expr$expression, c(1-quantiles,quantiles), na.rm = T)[2], "low",
                             ifelse(expr$expression >= quantile(expr$expression, c(1-quantiles,quantiles), na.rm = T)[1],
                                    "high", "none"))

        #q <- paste("anno", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        q2 <- paste("DE", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        q3 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        q4 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")

        anno_tmp <- merge(x= expr[expr$group != "none",], y= data[[i]][["anno"]][[j]], by="row.names")

        data_tmp <- data[[i]][["data"]][[j]][, anno_tmp$Row.names]

        sample <- as.factor(anno_tmp$group)
        design.mat <- model.matrix(~0+sample)
        colnames(design.mat) <- levels(sample)

        contrast.matrix <- makeContrasts(Diff= low - high, levels = design.mat)

        fit <- lmFit (data_tmp, design.mat)
        fit <- contrasts.fit(fit, contrast.matrix)
        fit <- eBayes(fit)

        DE_table <- topTable(fit, coef = "Diff", p.value = 1, adjust.method = adjust.method, lfc = log2(1), number = 100000)
        rank <- fit$coefficients[order(fit$coefficients[,1],decreasing = T),]

        gse <- suppressWarnings(fgsea(pathways = gs_list,
                                      stats = rank,
                                      minSize=1,
                                      maxSize=Inf,
                                      nperm=1000))

        container[[i]][["anno"]][[paste(i,j, sep = "_")]] <- anno_tmp
        container[[i]][["data"]][[paste(i,j, sep = "_")]] <- data_tmp
        container[[i]][["DE_genes"]][[paste(i,j, sep = "_")]] <- DE_table
        container[[i]][["Rank_genelist"]][[paste(i,j, sep = "_")]] <- rank
        container[[i]][["gsea"]][[paste(i,j, sep = "_")]] <- gse

        if (summ==T) {

          container[["summary"]][["Rank"]][[paste(i,j, sep = "_")]] <- rank
          container[["summary"]][["gsea"]][[paste(i,j, sep = "_")]] <- gse

          # Add the summary of the anno here
          tmp <- anno_tmp

          cont_anno <- list()

          for (n in colnames(tmp)[-c(1,2,3)]) {

            if (is.numeric(tmp[[n]])==T) {
              list <- tmp[,c("group", n)]
              list <- t.test(tmp[tmp$group=="high",][[n]], tmp[tmp$group=="low",][[n]], paired = F, var.equal = T, )
            }
            if (is.numeric(tmp[[n]])==F) {
              list <- tmp[,c("group", n)]
              list <- table(list)
              list <- prop.table(list, margin = 1)*100
            }

            cont_anno[[n]] <- list
          }

          container[["summary"]][["anno"]][[paste(i,j, sep = "_")]] <- cont_anno

        }

        if (length(names(data[[i]][["metadata"]]))>0) {

          for (k in names(data[[i]][["metadata"]])) {

            tmp_metadata <- merge(expr[expr$group != "none",], data[[i]][["metadata"]][[k]], by="row.names")

            container[[i]][["metadata"]][[paste(i, j, k, sep = "_")]] <- tmp_metadata

            if (summ==T) {

              tmp_metadata$expression <- NULL
              tmp_metadata <- suppressMessages(melt(tmp_metadata))
              tmp_metadata2 <- summarySE(tmp_metadata, groupvars = c("group", "variable"), measurevar = "value", na.rm = T)[,c(1,2,4)]
              tmp_metadata2 <- merge(tmp_metadata2[tmp_metadata2$group=="high",], tmp_metadata2[tmp_metadata2$group=="low",], by= "variable")
              tmp_metadata2 <- data.frame(variable=tmp_metadata2$variable, high_mean=tmp_metadata2$value.x, low_mean=tmp_metadata2$value.y, fc_low_high=tmp_metadata2$value.y/tmp_metadata2$value.x)

              # Calculate the pvalue
              tmp_metadata_pval <- compare_means(value~group, data = tmp_metadata, method = "t.test", paired = F, group.by = "variable", var.equal = FALSE, p.adjust.method = "none")

              # Merging with the sign value
              tmp_metadata <- merge(tmp_metadata2, tmp_metadata_pval[,c(1,5,8)], by = "variable")

              container[["summary"]][["metadata"]][[paste(i, j, k, sep = "_")]] <- tmp_metadata

            }

          }

        }

      }

      else {
        print(paste(gene, "is not present in", j, sep = " "))
      }

    }

  }

  class(container) <- "huva_experiment"

  return(container)
}

#' Summarise huva experiment in list or string
#'
#' @param huva_experiment huva_experiment class object calculate with the run_huva experiment function.
#' @param include_gene LOGIC. If TRUE genes fold change will be included in the summary
#' @param include_gsea LOGIC. If TRUE NES and p value from GSEA will be included in the analysis
#' @param include_metadata LOGIC. If TRUE metadata fold change in the two groups will be included in the summary
#' @param gene_list Filter for genes to be included in the summary (charachter vector)
#' @param metadata_list Filter for metadata tables to be included in the anaysis
#' @param one_line_result LOGIC. If TRUE the summary will be transformed in a single string
#' @return summary of huva_experiment
#' @examples
#'
#' 2+2
#'
#' @import limma
#' @import fgsea
#' @import Rmisc
#' @import ggpubr
#' @import reshape2
#' @export
summary_huva_experiment <- function(huva_experiment, include_gene=T, include_gsea=T, include_metadata=F, gene_list , metadata_list=NULL,one_line_result=F) {

  if (class(huva_experiment)!="huva_experiment") {
    print("use huva_experiment class object for reliable results")
  }

  container <- list()

  if ("summary" %in% names(huva_experiment)) {

    if (include_gene==T) {

      data <- data.frame(gene_name=gene_list)
      rownames(data) <- data$gene_name

      for (i in names(huva_experiment$summary$Rank)) {
        #print(i)
        data_2 <- data.frame(huva_experiment$summary$Rank[[i]])
        colnames(data_2)[1] <- i
        data_2$gene_name <- names(huva_experiment$summary$Rank[[i]])
        data <- merge(data, data_2, by = "gene_name", all.x = T)
      }

      rownames(data) <- data$gene_name
      data$gene_name <- NULL
      data$mean <- rowMeans(data, na.rm = T)

      container[["data"]] <- data

    }

    if (include_gsea==T) {

      gsea <- as.data.frame(huva_experiment$summary$gsea[[1]]$pathway)
      colnames(gsea) <- "pathway"

      for (i in 1:length(huva_experiment$summary$gsea)) {
        #print(names(huva_experiment$summary$gsea[i]))
        tmp3 <- huva_experiment$summary$gsea[[i]][, c("pathway", "NES")]
        gsea <- merge(gsea, tmp3, by="pathway")
        colnames(gsea)[i+1] <- names(huva_experiment$summary$gsea[i])
      }

      rownames(gsea) <- gsea$pathway
      gsea$pathway <- NULL
      gsea$mean <- rowMeans(gsea, na.rm = T)

      container[["gsea"]] <- gsea

    }

    if (include_metadata==T) {

      if (is.null(metadata_list)) {

        container[["metadata"]] <- huva_experiment[["summary"]][["metadata"]]

      }

      else {

        container[["metadata"]] <- huva_experiment[["summary"]][["metadata"]][metadata_list]

      }

    }

    if (one_line_result==T) {

      one_line <- vector()

      if (include_gene==T) {

        data_l <- data$mean
        names(data_l) <- rownames(data)

        one_line <- c(one_line, data_l)

      }

      if (include_gsea==T) {
        gsea_l <- gsea$mean
        names(gsea_l) <- rownames(gsea)

        one_line <- c(one_line, gsea_l)
      }

      if (include_metadata==T) {

        for (i in names(container[["metadata"]])) {

          meta_l <- container[["metadata"]][[i]][["fc_low_high"]]
          names(meta_l) <- container[["metadata"]][[i]][["variable"]]

          meta_l_stat <- container[["metadata"]][[i]][["p"]]
          names(meta_l_stat) <- paste("pval ",container[["metadata"]][[i]][["variable"]], sep = "")

          one_line <- c(one_line, meta_l, meta_l_stat)

        }
      }

      return(one_line)

    }

    else {

      class(container) <- "huva_experiment_summary"

      return(container)

    }



  }

  else {

    print("huva_experimet do not contain summarized results")

  }
}

#' Calculate huva experiment on the variance of a seleted metadata paramiter
#'
#' @param data huva_dataset class object.
#' @param phenotype gene name used for the analysis
#' @param quantiles definition of the quantile of segregation of the samples, quantiles are always simmetrical
#'     between high and low groups (e.g. \code{quantile 0.1} will use the 10th and 90th percentiles)
#' @param gs_list class list object defining gene sets to be included in the analysis (to generate this file see
#'     the documentation of fgse)
#' @param metadata_table name of the metadata table in wich the selected paramiter can be found
#' @param study defines the study to be used in the analysis
#' @param summ default TRUE, define if the summary of the huva experiment will be calculated
#' @param adjust.method p value adjustment method used to corret the DE genes analysis
#' @return huva_experiment based on the segregaation of the phenotype
#' @examples
#'
#' 2+2
#'
#' @import limma
#' @import fgsea
#' @import Rmisc
#' @import ggpubr
#' @import reshape2
#' @import stats
#' @export
run_huva_phenotype <- function(data, phenotype, study, metadata_table, quantiles, gs_list, summ=T, adjust.method="none") {

  if (class(data)!="huva_dataset") {
    print("Use huva_dataset class object for reliable results")
  }

  container <- list()

  print(paste("Binning on ", phenotype, sep = ""))

  print(paste("Using the ", study, " dataset for the analysis", sep = ""))

  data <- data[[study]] # selecton of only the study used for the analysis

  if (phenotype %in% colnames(data[["metadata"]][[metadata_table]])==F) {

    print("selected phenotype not found in the metadata table")

  }

  for (j in names(data[["data"]])) {

    binning <- as.data.frame(data[["metadata"]][[metadata_table]][,phenotype])
    rownames(binning) <- rownames(data[["metadata"]][[metadata_table]])
    colnames(binning) <- "phenotype"

    binning$group <- ifelse(binning$phenotype <= quantile(binning$phenotype, c(1-quantiles, quantiles), na.rm = T)[2], "low",
                            ifelse(binning$phenotype >= quantile(binning$phenotype, c(1-quantiles, quantiles), na.rm = T)[1],
                                   "high", "none"))

    #q <- paste("anno", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
    #q2 <- paste("DE", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
    #q3 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
    #q4 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")

    anno_tmp <- merge(x=binning[binning$group != "none",], y=data[["anno"]][[j]], by="row.names")

    data_tmp <- data[["data"]][[j]][, anno_tmp$Row.names]

    sample <- as.factor(anno_tmp$group)
    design.mat <- model.matrix(~0+sample)
    colnames(design.mat) <- levels(sample)

    contrast.matrix <- makeContrasts(Diff= low - high, levels = design.mat)

    fit <- lmFit (data_tmp, design.mat)
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)

    DE_table <- topTable(fit, coef = "Diff", p.value = 1, adjust.method = adjust.method, lfc = log2(1), number = 100000)
    rank <- fit$coefficients[order(fit$coefficients[,1],decreasing = T),]

    print(paste("GSEA on ", j, sep = ""))

    gse <- suppressWarnings(fgsea(pathways = gs_list,
                                  stats = rank,
                                  minSize=1,
                                  maxSize=Inf,
                                  nperm=1000))

    container[[study]][["anno"]][[paste(study, j, sep = "_")]] <- anno_tmp
    container[[study]][["data"]][[paste(study, j, sep = "_")]] <- data_tmp
    container[[study]][["DE_genes"]][[paste(study, j, sep = "_")]] <- DE_table
    container[[study]][["Rank_genelist"]][[paste(study, j, sep = "_")]] <- rank
    container[[study]][["gsea"]][[paste(study, j, sep = "_")]] <- gse

    if (summ==T) {

      container[["summary"]][["Rank"]][[paste(study, j, sep = "_")]] <- rank
      container[["summary"]][["gsea"]][[paste(study, j, sep = "_")]] <- gse

      # Add the summary of the anno here
      tmp <- anno_tmp

      cont_anno <- list()

      for (n in colnames(tmp)[-c(1,2,3)]) {

        if (is.numeric(tmp[[n]])==T) {
          list <- tmp[,c("group", n)]
          list <- t.test(tmp[tmp$group=="high",][[n]], tmp[tmp$group=="low",][[n]], paired = F, var.equal = T, )
        }
        if (is.numeric(tmp[[n]])==F) {
          list <- tmp[,c("group", n)]
          list <- table(list)
          list <- prop.table(list, margin = 1)*100
        }

        cont_anno[[n]] <- list
      }

      container[["summary"]][["anno"]][[paste(study, j, sep = "_")]] <- cont_anno

    }

  }

  if(length(names(data[["metadata"]]))>0) {

    print("analysis of metadata")

    for (k in names(data[["metadata"]])) {

      tmp_metadata <- merge(binning[binning$group != "none",], data[["metadata"]][[k]], by="row.names")

      container[[study]][["metadata"]][[paste(study, j, k, sep = "_")]] <- tmp_metadata

      if (summ==T) {

        tmp_metadata$expression <- NULL
        tmp_metadata <- suppressMessages(melt(tmp_metadata))
        tmp_metadata2 <- summarySE(tmp_metadata, groupvars = c("group", "variable"), measurevar = "value", na.rm = T)[,c(1,2,4)]
        tmp_metadata2 <- merge(tmp_metadata2[tmp_metadata2$group=="high",], tmp_metadata2[tmp_metadata2$group=="low",], by= "variable")
        tmp_metadata2 <- data.frame(variable=tmp_metadata2$variable, high_mean=tmp_metadata2$value.x, low_mean=tmp_metadata2$value.y, fc_low_high=tmp_metadata2$value.y/tmp_metadata2$value.x)

        # Calculate the pvalue
        tmp_metadata_pval <- compare_means(value~group, data = tmp_metadata, method = "t.test", paired = F, group.by = "variable", var.equal = FALSE, p.adjust.method = "none")

        # Merging with the sign value
        tmp_metadata <- merge(tmp_metadata2, tmp_metadata_pval[,c(1,5,8)], by = "variable")

        container[["summary"]][["metadata"]][[paste(study, j, k, sep = "_")]] <- tmp_metadata

      }


    }

  }

  class(container) <- "huva_experiment"

  return(container)


}

#' Calculate huva experiment on variance in the signature enrichment of a provided gene set
#'
#' @param data huva_dataset class object.
#' @param gene_set signature vector to be used for the definition of high and low groups
#' @param GSVA.method method for single sample signature enrichment, see ?gsva for options
#' @param quantiles definition of the quantile of segregation of the samples, quantiles are always simmetrical
#'     between high and low groups (e.g. \code{quantile 0.1} will use the 10th and 90th percentiles)
#' @param gs_list class list object defining gene sets to be included in the analysis (to generate this file see
#'     the documentation of fgse)
#' @param summ bolean. define if the summary of the huva experiment will be calculated
#' @param datasets_list characher vector used to filter the dataset in the data objet for the analysis
#' @param adjust.method p value adjustment method used to corret the DE genes analysis
#' @return huva_experiment
#' @examples
#'
#' 2+2
#'
#' @import limma
#' @import fgsea
#' @import Rmisc
#' @import ggpubr
#' @import reshape2
#' @import GSVA
#' @export
run_huva_signature <- function(data=huva_default_dataset, gene_set, quantiles, gs_list,summ=T, datasets_list=NULL, adjust.method="none", GSVA.method="gsva") {

  gs <- list()
  gs$bin_gs <- gene_set

  if (class(data)!="huva_dataset") {
    error("Use huva_dataset class object to run the huva_experiment function")
  }

  container <- list()

  print(paste("Binning on selected signature", sep = ""))

  if (is.null(datasets_list)==F) {
    # This is to allow the selection of the datasets to use in the analysis
    data <- data[datasets_list]
  }

  for (i in names(data)) {

    for (j in names(data[[i]][["data"]])) {

      expr <- suppressWarnings(gsva(data[[i]][["data"]][[j]], gset.idx.list = gs, method=GSVA.method))
      expr <- as.data.frame(t(expr))
      colnames(expr) <- c("expression")

      # Calculation of the percentiles

      expr$group <- ifelse(expr$expression<= quantile(expr$expression, c(1-quantiles,quantiles),na.rm = T)[2], "low",
                           ifelse(expr$expression >= quantile(expr$expression, c(1-quantiles,quantiles), na.rm = T)[1],
                                  "high", "none"))

      # q <- paste("anno", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
      # q2 <- paste("DE", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
      # q3 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
      # q4 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")

      anno_tmp <- merge(x= expr[expr$group != "none",], y= data[[i]][["anno"]][[j]], by="row.names")

      data_tmp <- data[[i]][["data"]][[j]][, anno_tmp$Row.names]

      sample <- as.factor(anno_tmp$group)
      design.mat <- model.matrix(~0+sample)
      colnames(design.mat) <- levels(sample)

      contrast.matrix <- makeContrasts(Diff= low - high, levels = design.mat)

      fit <- lmFit (data_tmp, design.mat)
      fit <- contrasts.fit(fit, contrast.matrix)
      fit <- eBayes(fit)

      DE_table <- topTable(fit, coef = "Diff", p.value = 1, adjust.method = adjust.method, lfc = log2(1), number = 100000)
      rank <- fit$coefficients[order(fit$coefficients[,1],decreasing = T),]

      gse <- suppressWarnings(fgsea(pathways = gs_list,
                                    stats = rank,
                                    minSize=1,
                                    maxSize=Inf,
                                    nperm=1000))

      container[[i]][["anno"]][[paste(i, j, sep = "_")]] <- anno_tmp
      container[[i]][["data"]][[paste(i, j, sep = "_")]] <- data_tmp
      container[[i]][["DE_genes"]][[paste(i, j, sep = "_")]] <- DE_table
      container[[i]][["Rank_genelist"]][[paste(i, j, sep = "_")]] <- rank
      container[[i]][["gsea"]][[paste(i, j, sep = "_")]] <- gse

      if (summ==T) {

        container[["summary"]][["Rank"]][[paste(i, j, sep = "_")]] <- rank
        container[["summary"]][["gsea"]][[paste(i, j, sep = "_")]] <- gse

        # Add the summary of the anno here
        tmp <- anno_tmp

        cont_anno <- list()

        for (n in colnames(tmp)[-c(1,2,3)]) {

          if (is.numeric(tmp[[n]])==T) {
            list <- tmp[,c("group", n)]
            list <- t.test(tmp[tmp$group=="high",][[n]], tmp[tmp$group=="low",][[n]], paired = F, var.equal = T, )
          }
          if (is.numeric(tmp[[n]])==F) {
            list <- tmp[,c("group", n)]
            list <- table(list)
            list <- prop.table(list, margin = 1)*100
          }

          cont_anno[[n]] <- list
        }

        container[["summary"]][["anno"]][[paste(i, j, sep = "_")]] <- cont_anno

      }

      if (length(names(data[[i]][["metadata"]]))>0) {

        for (k in names(data[[i]][["metadata"]])) {

          tmp_metadata <- merge(expr[expr$group != "none",], data[[i]][["metadata"]][[k]], by="row.names")

          container[[i]][["metadata"]][[paste(i, j, k, sep = "_")]] <- tmp_metadata

          if (summ==T) {

            tmp_metadata$expression <- NULL
            tmp_metadata <- suppressMessages(melt(tmp_metadata))
            tmp_metadata2 <- summarySE(tmp_metadata, groupvars = c("group", "variable"), measurevar = "value", na.rm = T)[,c(1,2,4)]
            tmp_metadata2 <- merge(tmp_metadata2[tmp_metadata2$group=="high",], tmp_metadata2[tmp_metadata2$group=="low",], by= "variable")
            tmp_metadata2 <- data.frame(variable=tmp_metadata2$variable, high_mean=tmp_metadata2$value.x, low_mean=tmp_metadata2$value.y, fc_low_high=tmp_metadata2$value.y/tmp_metadata2$value.x)

            # Calculate the pvalue
            tmp_metadata_pval <- compare_means(value~group, data = tmp_metadata, method = "t.test", paired = F, group.by = "variable", var.equal = FALSE, p.adjust.method = "none")

            # Merging with the sign value
            tmp_metadata <- merge(tmp_metadata2, tmp_metadata_pval[,c(1,5,8)], by = "variable")

            container[["summary"]][["metadata"]][[paste(i, j, k, sep = "_")]] <- tmp_metadata

          }

        }

      }

    }

  }

  class(container) <- "huva_experiment"

  return(container)
}
