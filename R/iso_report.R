
#' Gene multi-panel report, including gene locus and junction quantification,
#' average UMI counts per isoform (p1) and dotplot on normalized scaled data (p2).
#'
#' @param obj Seurat object
#' @param obj_assay Assay to pull data from
#' @param marker_list Output from ISO_SWITCH_ALL call
#' @param gene Gene to report on
#' @param gtf_df Data frame extract from gene annotation file
#' @param transcript_meta Transcript metadata
#'
#' @return ggplot object
#' @export
#' @importFrom magrittr %>%
#' @import Seurat
#' @import dplyr
#' @import ggplot2
#' @import forcats
#' @import patchwork
#' @import glue

#'
#' @examples
isoswitch_report <- function(obj, obj_assay, marker_list, gene, gtf_df, transcript_metadata) {

  normal_data <- obj@assays[[obj_assay]]@data
  isofs <- grep(paste("^", gene,"\\.", sep=""), rownames(normal_data), value=TRUE)
  assert_that(length(isofs) > 0, msg = "gene not found")

  # metadata shared by all panels (feature => short_name, color, order)
  meta <- ._build_plot_metadata(isofs, normal_data, transcript_metadata)

  # [Locus] _________________
  n_isofs = length(isofs)
  loc_plot <- ._isoswitch_report.locus(gtf_df, gene, meta)
  loc_h <- if(n_isofs <= 3) 1.3 else n_isofs * 0.34

  # [Junction] _________________

  if(!is.null(obj@assays$junction)) {
    jct_plot <- ._isoswitch_report.junctions(obj, "junction", gtf_df, gene, meta)
    jct_h <- 1.0
  } else {
    jct_plot <- plot_spacer()
    jct_h <- 0.01
  }

  # [P1 _________________ { umi counts }
  p1 <- ._isoswitch_report.umi_counts(obj, obj_assay, gene, meta)
  p1_celltype_order <- levels(p1$data$cell_type)

  # [P2]_________________ { dotplot }
  # use same cell_type order in y axis as in p1
  p2 <- ._isoswitch_report.dotpot(obj, obj_assay, meta, celltype_order=p1_celltype_order)

  # [PATCHWORK ] _________________
  pw <- loc_plot /
    jct_plot  /
    plot_spacer() /
    (p1 | p2)

  pw <- pw +
    plot_annotation(title = gene,
                    theme = theme(plot.title = element_text(size = 30))) +
    # 'null' expands to available space
    plot_layout(heights = unit(c(loc_h,  jct_h,  0.15,      5),
                               c( 'cm',   'cm',  'cm', 'null')))


  return(pw)

}
# ______________________________________________________________________________

#' Stripped-down version of isoswitch_report with only barplot
#'
#' @param obj Seurat object
#' @param obj_assay Assay to pull data from
#' @param marker_list Output from ISO_SWITCH_ALL call
#' @param gene Gene to report on
#' @param gtf_df Data frame extract from gene annotation file
#' @param transcript_meta Transcript metadata

#'
#' @return tbc
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @import forcats
#' @import patchwork
#' @import glue
#'
#' @examples
isoswitch_report_short <- function(obj, obj_assay, marker_list, gene, transcript_metadata) {

  normal_data <- obj@assays[[obj_assay]]@data
  isofs <- grep(paste("^", gene,"\\.", sep=""), rownames(normal_data), value=TRUE)
  assert_that(length(isofs) > 0, msg = "gene not found")

  # metadata shared by all panels (feature => short_name, color, order)
  meta <- ._build_plot_metadata(isofs, normal_data, transcript_metadata)

  # [P1 _________________ { umi counts }
  p1 <- ._isoswitch_report.umi_counts(obj, obj_assay, gene, meta, legend=TRUE)

  # [PATCHWORK ] _________________
  # extract the legend
  p1_legend <- cowplot::get_legend(p1)
  # plot p1 and legend together
  pw <- ((p1 + theme(legend.position = 'none')) / p1_legend) +
    plot_layout(heights = c(10, 1))

  return(pw)
}

# ______________________________________________________________________________

#' Returns a dataframe of metadata shared by all panels (feature => short_name, color, order)
#'
#' @param isofs TBC
#' @param normal_data TBC
#' @param transcript_meta TBC
#'
#' @return Returns a dataframe of metadata shared by all panels (feature => short_name, color, order)
#' @export
#'
#' @examples
#'
._build_plot_metadata <- function(isofs, normal_data, transcript_meta) {

  # order by expression
  ordered_isofs <- data.frame(isofs = isofs, expr = Matrix::rowSums(normal_data[isofs, ])) %>%
    arrange(desc(expr)) %>%
    pull(isofs)

  # custom scale, brewer.pal(n=10, name="Set3") reordered
  custom_colors = c("#FB8072","#80B1D3","#8DD3C7","#FFFFB3","#BEBADA",
                    "#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD")

  # build feature metadata (feature => short_name, color, order) shared by all panels
  meta <- data.frame( feature = ordered_isofs ) %>%
    separate(feature, into=c("gene_id", "transcript_id"), sep="\\.\\.", remove=FALSE) %>%
    left_join(transcript_meta, by=c("transcript_id"="ensembl_transcript_id")) %>%
    mutate(color = custom_colors[1:length(isofs)])

  return(meta)
}



# ______________________________________________________________________________




#' Title
#'
#' @param obj tbc
#' @param obj_assay tbc
#' @param meta tbc
#' @param celltype_order tbc
#' @param switch tbc
#'
#' @return tbc
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @import forcats
#' @import patchwork
#' @import glue
#'
#' @examples
._isoswitch_report.dotpot <- function(obj, obj_assay, meta, celltype_order=NULL, switch=NULL) {

  normal_data <- obj@assays[[obj_assay]]@data
  scaled_data <- obj@assays[[obj_assay]]@data

  df <- data.frame()
  for(cell_type in levels(obj@active.ident)) {


    cell_list <- Seurat::WhichCells(obj, ident = cell_type)
    n_cells = length(cell_list)
    isofs <- meta$feature

    normal_counts <- normal_data[isofs, cell_list, drop=FALSE]

    cell_df <- data.frame(feature = rownames(normal_counts)) %>%
      separate(feature, into=c("gene_id", "transcript_id"), sep="\\.\\.", remove=FALSE) %>%
      mutate(avg_scaled_value = Matrix::rowSums(scaled_data[isofs, cell_list, drop=FALSE]) / n_cells,
             cell_expr  = Matrix::rowSums(normal_counts > 0),
             cell_count = n_cells,
             type = cell_type)

    df <- rbind(df, cell_df)
  }

  # Apply same logic on which cell_types to show as ._isoswitch_report.umi_counts
  # -> cell_types with NO counts are removed
  # -> For cts with counts < 1%, keep them but set to 0 to avoid drawing dot on dotplot
  d2 <- df %>%
    filter(cell_expr > 0 ) %>%
    mutate(perc_expr = (cell_expr/cell_count)*100) %>%
    mutate(perc_expr = if_else(perc_expr < 1, 0, perc_expr)) %>%
    left_join(meta, by="transcript_id")

  # limits for the expression scale, ensures it's symmetric around zero
  top_value <- max(abs(df$avg_scaled_value))

  # main dotplot figure
  p2 <- ggplot(d2,
               aes(x=factor(external_transcript_name, levels=meta$external_transcript_name),
                   y=factor(type, levels=celltype_order),
                   fill = avg_scaled_value,
                   size = perc_expr)) +
    geom_point(shape=21, color="#333333") +
    scale_fill_distiller(palette="RdBu", limits=c(0, top_value)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=0, hjust=0.5, size=8),
          axis.text.y = element_text(size=8),
          legend.title = element_text(size=9),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank()) +
    labs(x=NULL, y=NULL, fill="Expression", size="% cells")

  # add switch segment
  # p2 <- p2 + geom_segment(x=meta[meta$transcript_id==switch$t1, "external_transcript_name"], y=switch$c1,
  #                         xend=meta[meta$transcript_id==switch$t2, "external_transcript_name"], yend=switch$c2,
  #                         alpha=1, linetype=5, color="gray", lineend="round", size=0.1)
  return(p2)
}

# ______________________________________________________________________________
#' Title
#'
#' @param gtf_df tbc
#' @param gene tbc
#' @param meta tbc
#' @param legend tbc
#'
#' @return tbc
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @import ggtranscript
#'
#' @examples
._isoswitch_report.locus <- function(gtf_df, gene, meta, legend=TRUE){

  # get exons for this gene (only transcripts shown)
  exons <- gtf_df %>%
    filter(gene_name==gene,
           transcript_id %in% meta$transcript_id,
           type=="exon") %>%
    left_join(meta, by="transcript_id")

  # force color mapping for each feature
  manual_colors <- tibble::deframe(select(meta, external_transcript_name, color))

  p1 <- ggplot(exons, aes(xstart = start, xend = end, y = factor(external_transcript_name, levels=meta$external_transcript_name))) +
    geom_range(aes(fill = external_transcript_name) ) +
    scale_fill_manual(values = manual_colors) +
    scale_y_discrete(limits=rev) +
    geom_intron(data = to_intron(exons, "transcript_id"), aes(strand = strand)) +
    labs(y=NULL, fill="Transcript") +
    theme_minimal()

  if(legend)
    p1 <- p1 +  theme(legend.key.height= unit(0.15, 'cm'),
                      legend.key.width= unit(0.40, 'cm'))
  else
    p1 <- p1 + theme(legend.position = "none")

  return(p1)
}

# ______________________________________________________________________________

#' Title
#'
#' @param obj tbc
#' @param assay tbc
#' @param gtf_df tbc
#' @param gene tbc
#' @param meta tbc
#'
#' @return tbc
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @import ggtranscript
#'
#' @examples
._isoswitch_report.junctions <- function(obj, obj_assay, gtf_df, gene, meta){

  junction_counts <- obj@assays[[obj_assay]]@counts
  idx <- grep(paste("^", gene,"\\.", sep=""), rownames(junction_counts), value=TRUE)
  gene_junctions <- junction_counts[idx, ]

  # chr/strand lookup
  gene_location <- gtf_df %>%
    dplyr::select(gene_name, seqnames, strand) %>%
    filter(gene_name==gene) %>% head(1)

  # get junction data (start, end, quantif), keep only jcts over 0.1 UMI/cell
  jct_df <- data.frame(junction = rownames(gene_junctions)) %>%
    separate(junction, into=c("gene_id","coords"), sep="\\.\\.") %>%
    separate(coords, into=c("start","end"), sep="-", convert=TRUE) %>%
    mutate(avg = Matrix::rowSums(gene_junctions)/ncol(gene_junctions)) %>%
    left_join(gene_location, by=c("gene_id"="gene_name")) %>%
    mutate(transcript_id = "ALL_EXONS")%>%
    filter(avg >= 0.1)

  # collapse exons for this gene (only transcripts shown) into a fake single id
  exons <- gtf_df %>%
    filter(gene_name==gene,
           type=="exon",
           transcript_id %in% meta$transcript_id) %>%
    mutate(transcript_id = "ALL_EXONS")

  ggplot(exons, aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range() +
    geom_intron(data = to_intron(exons, "transcript_name"), aes(strand = strand)) +
    geom_junction(data = jct_df, junction.y.max = 0.55) +
    geom_junction_label_repel(data = jct_df, aes(label = round(avg, 1)), force=3, junction.y.max = 1, label.size=0.12, size=2.4) +
    theme_minimal() +
    labs(y=NULL)
}


# ______________________________________________________________________________

#' Stacked barplot representing average raw transcript counts per cell for each cell type
#'
#' @param obj Seurat object
#' @param gene Gene
#' @param meta Metadata computed by isoswitch_report function
#' @param legend Legend boolead, defaults to False for iso_dotplo2 reports
#'
#' @return ggplot object
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @import forcats
#' @import patchwork
#' @import glue
#'
#' @examples
._isoswitch_report.umi_counts <- function(obj, obj_assay, gene, meta, legend=FALSE) {

  # get celltype_features data for each cell type
  df <- map_df(levels(obj@active.ident),
               celltype_features, obj=obj, count_matrix=obj@assays[[obj_assay]]@counts, features=meta$feature)

  # sort cell types by expresion of major isoform - reuse sorting in meta
  celltype_order <- filter(df, transcript_id == meta$transcript_id[[1]]) %>%
    arrange(avg) %>%
    pull(cell_type)

  # relevel factors, isoforms by average total UMI sum
  df <- mutate(df,
               transcript_id = forcats::fct_reorder(transcript_id, as.numeric(avg), .fun = mean),
               cell_type = fct_relevel(cell_type, levels = celltype_order))

  # force color mapping for each feature
  color_mapping <- select(meta, transcript_id, color) %>% tibble::deframe()

  p1<- ggplot(df, aes(y=cell_type, x=avg, fill=transcript_id)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = color_mapping, labels=meta$external_transcript_name) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          axis.title.x = element_text(size=8)) +
    guides(fill = guide_legend(nrow = 1, reverse = TRUE)) +
    labs(x="Avg UMI per cell", y=NULL, fill="")

  if(legend)
    p1 <- p1 + theme(legend.position = "bottom",
                     legend.text = element_text(size=7),
                     legend.key.size = unit(0.6,"line"),
                     legend.margin=margin(0,0,0,0),
                     legend.box.margin=margin(0,0,0,0))
  else
    p1 <- p1 + theme(legend.position = "none")

  return(p1)
}
