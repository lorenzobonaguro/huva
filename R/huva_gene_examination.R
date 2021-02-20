#' Examine the profile of expression of a selected gene and investigate the relationship with annotation parameters of
#' metadata
#'
#' @title Exploratory huva experiment
#' @description This funtion looks into the expression of a gene of interest (GOI) in the specified dataset or across other user-provided datasets.  
#' @param huva_dataset huva_dataset class object.
#' @param gene gene name used for the analysis.
#' @param bins graphical option, see ggplot2 for details.
#' @param alpha graphical option, see ggplot2 for details.
#' @return huva_gene_examination class object including a graphical overview of the distribution of gene expression of a selected GOI.
#' @seealso run_huva_experiment
#' @import ggplot2
#' @import reshape2
#' @examples
#' library(huva)
#' library(huva.db)
#' 
#' gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")
#'
#' @export
gene_exam <- function(huva_dataset = huva_default_data, gene = gene_name, bins = 30, alpha = 1){

  if (class(huva_dataset)!= "huva_dataset") {
    warning("Warning: For reliable results use class huva_dataset object")
  }

  print(paste("Analysis of", gene, "expression", sep = " "))

  container <- list()

  # If using all the datasets provided there is a loop to get the expression data of a selected gene and prepare table with all the annnotation and metadata for the analysis, all the
  # results will be stored in a list serving then as input for other functions

  for (i in names(huva_dataset)) {

    for (k in names(huva_dataset[[i]][["data"]])) {

      if (gene %in% rownames(huva_dataset[[i]][["data"]][[k]])) {

        dataframe <- huva_dataset[[i]][["data"]][[k]][gene,]
        dataframe <- melt(dataframe, value.name = "expression")

        container[["plot"]][[paste(i,k, sep = "_")]] <- ggplot(dataframe, aes(x = expression)) +
          geom_histogram(aes(y =..density..),bins = bins, alpha=alpha, fill="white", color="black") +
          geom_density(alpha =.2, color="black")+
          ggtitle(paste(gene, "-",k))+
          theme_minimal() +
          theme(aspect.ratio = 1)

        container[["expression"]][[paste(i,k, sep = "_")]] <- dataframe

        # Merging the expression with metadata

        container[["anno"]][[paste(i,k, sep = "_")]] <- merge(dataframe, huva_dataset[[i]][["anno"]][[k]], by="row.names",all.x=TRUE)

        if (length(names(huva_dataset[[i]][["metadata"]]))>0) {
          #print("Analysis of metadata for this dataset")

          for (l in names(huva_dataset[[i]][["metadata"]])) {

            container[["metadata"]][[paste(i, k, l, sep = "_")]] <- merge(container[["expression"]][[paste(i,k, sep = "_")]], huva_dataset[[i]][["metadata"]][[l]], by="row.names", all.x=TRUE)

          }

        }

      }

      else {
        print(paste(gene, "is not present in", k, sep = " "))
      }

    }
  }

  # Generating a combined expression table

  maxrow <- max(sapply(container$expression, nrow))

  mylist_mod <- lapply(container$expression, function(x,nRow){
    if(nrow(x) <  nRow){
      x[(nrow(x)+1):nRow,] <- NA
    }
    x
  }, nRow = maxrow)

  expression_sum <- as.data.frame(mylist_mod)

  colnames(expression_sum) <- names(container$expression)
  rownames(expression_sum) <- NULL

  container[["expression"]][["expression_summary"]] <- expression_sum

  container[["gene_name"]] <- gene

  # output
  
   class(container) <- "huva_gene_examination"

   return(container)

}
