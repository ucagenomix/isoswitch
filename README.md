
<!-- README.md is generated from README.Rmd. Please edit that file -->

# isoswitch

<!-- badges: start -->

<!-- badges: end -->

The goal of isoswitch is to facilitate the characterization of isoform
expression in long-read single-cell datasets (ScNaUmi-seq, Lebrigand et
al 2020). It includes a set of functions built on top of Seurat and
ggplot that can be used to search, visualize and document isoform switch
patterns.

## Installation

You can install the development version of isoswitch from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("atienza-ipmc/isoswitch")
```

## General Workflow

1.  Input data / object setup
2.  Isoform characterization
3.  Switch search
4.  Gene reports

### 1 - Input data / object setup

Isoswitch works with Seurat objects with gene- and isoform-level counts.

  - A standard gene-level “RNA” assay, with \[gene x cell\] matrix count
  - An isoform-level with \[isoform x cell\] matrix counts, By
    convention, the row names of the isoform count matrix follow the
    format “<gene_name>..<transcript_id> -\>”BCS1L..ENST00000359273"

### 2\. Isoform characterization

You can get a plot summarizing number of genes, number of transcripts,
distribution of isoforms and number of genes per cell type using the
method *plot\_assay\_stats*

``` r
plot_assay_stats(seurat, "isoform")
```

![alt text](./man/figures/Fig4_isosummary.png)

Counts

that containsic example which shows you how to solve a common problem:

``` r
library(isoswitch)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.svg" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
