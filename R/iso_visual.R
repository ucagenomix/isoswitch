

#
#' Plot summary stats on gene and isoform expression
#'
#' Produces a 4 panel plot including a) Number of total genes and transcripts
#' and percentage of genes expressing one or more isoforms. (b) Distribution
#' of the number of isoforms per gene. (c) Relationship between the number of
#' major isoforms (y-axis) and the total number of isoforms per gene (x-axis).
#' (d) Pseudo-bulk estimation of number of genes expressed per cell type
#' (each data point represents the aggregation of 10 cells)
#'
#'
#' @param obj Seurat object
#' @param assay Assay to pull data from
#' @param bin_size Bin size for pseudobulk aggregation of number of genes per cell
#'
#' @return ggplot object, patchwork of 4 plots
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @import forcats
#' @import patchwork
#'
#' @examples
plot_assay_stats <- function(obj, assay, bin_size=10) {

  iso_stats <- iso_compute_stats(obj@assays[[assay]]@counts, count_zeros=FALSE)

  # [p1 - Genes/Transcripts -------]
  total_genes = n_distinct(iso_stats$gene_id)
  singles = nrow(count(iso_stats, gene_id) %>% filter(n==1))
  single_label <- if(singles) paste0(round(100*singles/total_genes, 1), "%") else ""
  multis = total_genes - singles
  multis_label <- if(multis) paste0(round(100*multis/total_genes, 1), "%") else ""

  d1 <- data.frame()
  d1 <- rbind(d1, c("Genes", "mono-iso", singles, single_label))
  d1 <- rbind(d1, c("Genes", "multi-iso", multis, multis_label))
  d1 <- rbind(d1, c("Transcripts", "transcripts", n_distinct(iso_stats$transcript_id), ""))
  colnames(d1) <- c("t", "i", "n", "l")
  d1 <- mutate(d1, n = as.numeric(n),
               i = factor(i, levels = c("multi-iso", "mono-iso")))

  p1 <- ggplot(d1, aes(x=t, y=n, fill=i)) +
    geom_col() +
    geom_text(aes(label=l), colour="#444444", size=3, position = position_stack(vjust = .5)) +
    scale_fill_manual(values = c("#fcc729", "#337def", "#888888"), limits = c('mono-iso', 'multi-iso')) +
    theme_minimal() +
    theme(legend.position="bottom",
          legend.key.height = unit(4, "pt"),
          legend.key.width = unit(12, "pt")) +
    labs(x=NULL, y=NULL, fill="")

  # [p2 - Genes expressed per cell type  -------]
  celltype_data <- celltype_features_bulk(obj, assay, bulk_bin_size=bin_size)

  p2 <- ggplot(celltype_data, aes(x=reorder(t,n_genes), y=n_genes)) +
    geom_boxplot() +
    geom_jitter(size=0.5, alpha=0.4) +
    theme_minimal() +
    scale_y_continuous(limits=c(0, NA)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7.5),
          axis.title = element_text(size=8)) +
    labs(x = "Genes expressed per cell type",
         y= NULL)

  # [p3 - Nb isoforms per gene -------]

  bin_levels = c("1","2","3","4","5","6","7","8","9",">10")

  d3 <- count(iso_stats, gene_id) %>%
    count(n, name="nn") %>%
    mutate(bin = factor(if_else(n < 10, as.character(n), ">10"))) %>%
    mutate(bin = fct_relevel(bin, bin_levels)) %>%
    group_by(bin) %>%
    summarise(total=sum(nn)) %>%
    mutate(perc = 100 * total / total_genes)

  p3 <- ggplot(d3, aes(x=bin, y=total)) +
    geom_col() +
    geom_text(aes(label=glue("{round(perc,1)}%")), colour="#444444", vjust=-0.1, size=2.5) +
    scale_x_discrete(limits = bin_levels) +
    theme_minimal() +
    theme(axis.title = element_text(size=8)) +
    # leave room for the annotation on the top bar
    coord_cartesian(ylim = c(0, max(d3$total)),  clip = 'off') +
    xlab("Number of isoforms per gene") +
    ylab(NULL)

  # [p4 - Major isoform distribution boxplot -------]

  d4 <- iso_stats %>%
    group_by(n_isofs, gene_id) %>%
    summarise(n_majors = sum(is_major))

  p4 <- ggplot(d4, aes(x=n_isofs, y=n_majors)) +
    geom_point(alpha=0.1, size=0.25) +
    geom_jitter(size=0.25)+
    geom_smooth(method='lm', size=0.25) +
    # coord_cartesian(xlim=c(0,20), ylim=c(0,10)) +
    theme_minimal() +
    theme(axis.title = element_text(size=8),
          legend.title=element_text(size=8)) +
    labs(x = "Number of isoforms",
         y = "Major isoforms")

  pw <- (p1 / p3 / p4) | p2
  return(pw)
}



# ______________________________________________________________________________

#' Volcano-like plot showing log_f2c vs p-value for genes found in marker_list
#' returned by ISO_SWITCH_ALL
#'
#' @param obj Seurat object
#' @param markers List of markers returned by ISO_SWITCH_ALL
#' @param facet if TRUE, facet by individual cell type, defaults to FALSE
#' @param ncol specifies number of columns if facet=TRUE, defaults to 4
#' @param overlaps max number of overlapping labels, defaults to 10
#'
#' @return ggplot object
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
plot_marker_score <- function(obj, markers, facet=FALSE, ncol=4, overlaps=10){
  assert_that(nrow(markers) > 0, msg="Empty marker list")

  # keep only smallest p-value for each gene
  d1 <- if(facet) group_by(markers, cluster, geneId) else group_by(markers, geneId)
  d1 <- slice_min(d1, order_by=p_val_adj, n=1, with_ties=FALSE)

  # p-values=0 break the scaling, keep a very low p instead
  d1 <- mutate(d1, p_val = if_else(p_val==0, 1e-230, p_val))
  # rescale Y (p_vals) as Xs, and calculate distance to (0,0) - for color fill
  d1$y_rescaled = scales::rescale(-log10(d1$p_val), c(0, max(d1$avg_log2FC)))
  d1 <- mutate(d1, dist = sqrt(sum((c(avg_log2FC,y_rescaled)-c(0,0))^2)))

  p1 <- ggplot(d1, aes(x=avg_log2FC, y=p_val, color=dist)) +
    geom_point(alpha=1) +
    scale_y_continuous(trans=._reverselog_trans(10))  +
    scale_color_distiller(palette = "Reds", direction=1) +
    scale_x_continuous(breaks = seq(ceiling(min(markers$avg_log2FC)),floor(max(markers$avg_log2FC)),1))  +
    geom_text_repel(aes(label=geneId), size=2.6, max.overlaps=overlaps, color="#444444") +
    coord_cartesian(clip = 'off') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "#AAAAAA", size=0.7),
          axis.text.x = element_text(size=8),
          axis.title = element_text(size=8),
          legend.position = "none")

  if(facet){
    p1 <- p1 + facet_wrap(~cluster, ncol=ncol)
  }

  return(p1)
}

# ______________________________________________________________________________
#' Plot matrix of genes found per contrast
#'
#' @param obj Seurat object
#' @param marker_list Marker list returned by ISO_SWITCH_ALL
#'
#' @return ggplot object
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @export
#'
#' @examples
plot_marker_matrix <- function(obj, marker_list) {

  marker_list <-  group_by(marker_list, contrast) %>%
    summarise(n=n_distinct(geneId)) %>%
    tidyr::separate(contrast, into=c("c1", "c2"), sep=" vs ")

  p1<-ggplot(marker_list,
             aes(x=factor(c1, levels=levels(obj@active.ident)), y=factor(c2, levels=levels(obj@active.ident)), fill=n)) +
    geom_tile() +
    geom_text(aes(label=n), size=2) +
    theme_minimal() +
    scale_fill_distiller(palette="YlOrRd", direction=1) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          title = element_text(size=10),
          plot.title.position = "plot",
          legend.position='none') +
    labs(x=NULL, y=NULL, title = "Number of isoform-switching genes per contrast")

  return(p1)
}


# ______________________________________________________________________________
#' Plot matrix of genes found per contrast with row annotations
#'
#' @param obj Seurat object
#' @param marker_list Marker list returned by ISO_SWITCH_ALL
#'
#' @return ggplot object
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @export
#'
#' @examples
plot_marker_matrix_annotated <- function(obj, marker_list) {

  pairwise <-  group_by(marker_list, contrast) %>%
    summarise(n=n_distinct(geneId)) %>%
    separate(contrast, into=c("c1", "c2"), sep=" vs ")

  # separate annotations per row/column
  col_annot <- marker_list %>%
    separate(contrast, into=c("c1", "c2"), sep=" vs ")  %>%
    group_by(c1) %>%
    summarise(count_c=n_distinct(geneId))

  row_annot <- marker_list %>%
    separate(contrast, into=c("c1", "c2"), sep=" vs ")  %>%
    group_by(c2) %>%
    summarise(count_r=n_distinct(geneId))

  # zero-fill missing counts
  celltype_order <- levels(obj@active.ident)
  missing_cols <- data.frame(c1=setdiff(celltype_order, col_annot$c1), count_c=0)
  missing_rows <- data.frame(c2=setdiff(celltype_order, row_annot$c2), count_r=0)
  row_annot <- bind_rows(row_annot, missing_rows)
  col_annot <- bind_rows(col_annot, missing_cols)

  # combine (sum) the two
  combined_annot <- full_join(row_annot, col_annot, by=c("c2"="c1")) %>%
    mutate(total = count_r + count_c)

  pl <- ggplot(pairwise,
               aes(x=factor(c1, levels=celltype_order), y=factor(c2, levels=celltype_order))) +
    geom_tile(aes(fill=n)) +
    geom_text(aes(label=n), size=2) +
    theme_minimal() +
    scale_fill_distiller(palette="YlOrRd", direction=1) +
    scale_y_discrete(drop=FALSE)  +
    # scale_x_discrete(drop=FALSE)  +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          title = element_text(size=10),
          plot.title.position = "plot") +
    labs(x=NULL, y=NULL, title = "   ") +
    NoLegend()

  # add col/row annotation
  pl <- pl +
    #geom_text(data=col_annot, aes(x=c1, label=count_c), y=tail(layer_scales(pl)$y$range$range, n=1), vjust = -3.2, size = 2.5) +
    #geom_text(data=row_annot, aes(y=c2, label=str_pad(count_r,2)), x=tail(layer_scales(pl)$x$range$range, n=1), hjust = -2.5, size = 2.5) +
    geom_text(data=combined_annot, aes(y=c2, label=total), x=tail(layer_scales(pl)$x$range$range, n=1), hjust = -2.5, size = 2.5) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = unit(c(0.2,1,0.2,0.2), "cm"))

  # add total number of genes
  nb_unique_genes = n_distinct(marker_list$geneId)
  pl <- pl +
    annotation_custom(grid::textGrob(label = paste0("n=", nb_unique_genes),
                                     x = unit(0.90, "npc"),
                                     y = unit(0.09, "npc"),
                                     gp = grid::gpar(fontsize = 9.5)))
  return(pl)
}

# ______________________________________________________________________________

# -log10 scale -> volcano plots
._reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
            scales::log_breaks(base = base),
            domain = c(1e-100, Inf))
}
