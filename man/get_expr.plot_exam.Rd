% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/huva_explore_examination.R
\name{get_expr.plot_exam}
\alias{get_expr.plot_exam}
\title{get_expr.plot_exam}
\usage{
get_expr.plot_exam(
  huva_expression,
  study = "ALL",
  dataset = "ALL",
  bins = 30,
  alpha = 1
)
}
\arguments{
\item{huva_expression}{huva_expression class object}

\item{study}{character vector defining the studies to be used in the analysis (default = "ALL").}

\item{dataset}{character vector defining the datasets to be used in the analysis (default = "ALL").}

\item{bins}{graphical option, see ggplot2 for details (default = 30).}

\item{alpha}{graphical option, see ggplot2 for details (default = 1).}
}
\value{
ggplot class object, histogram reporting the expression of a selected GOI across selected datasets.
}
\description{
The function plots expression of the gene of interest (GOI) across all datasets to inspect possible
             distribution differences in a huva gene examination experiment.
}
\details{
Visualization of the expression histogram for a selected gene of interest
}
\examples{
library(huva)
library(huva.db)

gene_overview <- gene_exam(huva_dataset = huva.db, gene = "MYD88")

expr_exam.plot <- get_expr.plot_exam(huva_expression = gene_overview)

}
\seealso{
gene_exam, run_huva_experiment
}
