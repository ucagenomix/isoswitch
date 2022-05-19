
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

Below is a short overview of the package functionality using a fake
dataset.

### 1 - Input data / object setup

Isoswitch works with Seurat objects with gene- and isoform-level counts.

  - A standard gene-level “RNA” assay, with \[gene x cell\] matrix count
  - An isoform-level with \[isoform x cell\] matrix counts, By
    convention, the row names of the isoform count matrix follow the
    format “<gene_name>..<transcript_id> -\>”BCS1L..ENST00000359273"

### 2\. Isoform characterization

The method **iso\_compute\_stats** parses the isoform raw count matrix
and returns a data frame with the following structure:

``` r
stats <- iso_compute_stats(seurat@assays$multi@counts) %>% arrange(gene_id)
head(stats, n=4)
#>                 feature gene_id   transcript_id sum total_gene n_isofs max_sum
#> 1 A1BG..ENST00000596924    A1BG ENST00000596924   5          8       2       5
#> 2 A1BG..ENST00000598345    A1BG ENST00000598345   3          8       2       5
#> 3  A2M..ENST00000495709     A2M ENST00000495709  10         14       2      10
#> 4  A2M..ENST00000318602     A2M ENST00000318602   4         14       2      10
#>       perc is_major is_top
#> 1 62.50000     TRUE   TRUE
#> 2 37.50000     TRUE  FALSE
#> 3 71.42857     TRUE   TRUE
#> 4 28.57143    FALSE  FALSE
```

The method **plot\_assay\_stats** builds on this data to plot a summary
with number of genes, number of transcripts, distribution of isoforms
and number of genes per cell type.

``` r
plot_assay_stats(seurat, "isoform")
```

![alt text](./man/figures/Fig4_isosummary.png)

### 3\. Isoform switch search

The term “isoform switch” refers to an event where two isoforms of the
same gene are considered markers of different clusters.

The marker search is implemented on **ISO\_SWITCH\_ALL** which relies on
Seurat’s FindMarkers functionality.

``` r
clusters <- levels(seurat@active.ident)
switch_markers <- ISO_SWITCH_ALL(seurat, clusters, assay="isoform", min.pct=0, logfc.threshold=0.40, verbose=TRUE)
```

The result of **ISO\_SWITCH\_ALL** is a data frame of transcripts
considered statistically significant markers of certain clusters. The
method **compute\_switches** takes a list of markers as input and
combines them into switches, ranking the list of switches according to
different possible criteria.

``` r
switches <- compute_switches(switch_markers)
```

The helper method **format\_switch\_table** prettifies the switch table
for reports:

To facilitate the interpretation of this data, two methods
**plot\_marker\_matrix** and **plot\_marker\_score** generate plots with
1) heatmap of number of unique genes per contrast between clusters and
2) volcano-like showing p-values and average logFC for each gene with an
isoform switch

``` r
pl1 <- plot_marker_matrix(seurat, switch_markers) 
pl2 <- plot_marker_score(adult, switch_markers, facet=FALSE, overlaps=16)
pl1 | pl2 
```

![alt text](./man/figures/Fig7_isoswitch.png) Alternatively,
**plot\_marker\_score** can also plot individual plots for each cluster
analyzed

``` r
plot_marker_score(adult, switch_markers, facet=TRUE, ncol=3)
```

![alt text](./man/figures/Fig7_facet.png)

### 4\. Gene reports

After identifying genes of interest, **isoswitch\_report** produces a
detailed report of the isoform switch

``` r
isoswitch_report(seurat, "isoform", gene="HYAL2", marker_list=switch_markers) 
```

![alt text](./man/figures/Fig7_hyal2.png)
