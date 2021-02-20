# huva (R)
<img src="./logo/huva_logo_V3.png" width="20%" style="float: left;">

In the last decades, the spread of increasingly affordable transcriptomic techniques and the concomitant production of population-scale multi-layered datasets assembled extensive experimental data of different types on single healthy individuals. Within this population, most genetic variation and environmental factors are influencing gene expression with no clinical evidence of pathological states. Variance of gene expression, which is a characteristic of any given natural population, can be exploited to set up a conditional quasi loss- and gain-of-function “in population” experiment if expression values for the gene of interest (GOI) are available in a sufficiently large cohort. We describe here a novel tool, called *huva* (human variation), which takes advantage of population-scale multi-layer data to infer gene function and relationship between phenotype and gene expression. Within a reference dataset, *huva* derives two experimental groups, i.e. individuals with LOW or HIGH expression of the GOI, enabling the subsequent comparison of their transcriptional profile and functional parameters. We demonstrate that this approach robustly and efficiently identifies the phenotypic relevance of a GOI, allowing the stratification of genes according to shared functions. Additionally, we describe how *huva* can predict the phenotype of naturally occurring gain-of-function mutations in humans, further validating the biological power of the approach.

## how to install *huva*
You can istall *huva* from R with the following code:

```R
url <- "https_url_of_the_package"

devtools::install_git(url = url)
```

## Usage
For detailed informations usage check the *huva vignette* and the documentation of each function

## How to cite *huva*
If you use *huva* in your research project consider citing us: [huva: human variation in gene expression as surrogate for gene function; Bonaguro at al. 2021](weblink).

## Contact or follow us
For any problem of question regrding the *huva* package or this repositoy or you just want to be up to date on what is coming next, send us an [email](mailto:helphuva@uni-bonn.de) or follow us:  

<img src="./logo/twitter.png" width="12%" style="float: left;">  

[@LorenzoBonaguro](https://twitter.com/LorenzoBonaguro)  
[@AschenbrennerAC](https://twitter.com/AschenbrennerAC)  
[@LabSchultze](https://twitter.com/LabSchultze)