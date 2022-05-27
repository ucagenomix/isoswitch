

#' Find transcript markers of isoform switches
#'
#' @param seurat Seurat object
#' @param clusters Clusters to compare
#' @param assay Assay to pull data from
#' @param ... extra parameters passed on to Seurat FindMarkers call
#'
#' @return Marker list data frame
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import glue
#' @importFrom tibble rownames_to_column
#' @importFrom utils combn
#' @export
#'
#' @examples
ISO_SWITCH_ALL <- function(seurat, clusters, assay,  ... ){

  Seurat::DefaultAssay(seurat) <- assay
  get_markers <- function(cluster_i, cluster_j, obj)  {

    markers <- Seurat::FindMarkers(object = seurat, ident.1=cluster_i, ident.2=cluster_j, ...) %>%
      rownames_to_column("feature") %>%
      separate(feature, into=c("geneId", "transcriptId"), sep="\\.\\.", remove=FALSE) %>%
      mutate(cluster  = ifelse(avg_log2FC > 0, cluster_i, cluster_j),
             contrast = glue::glue("{cluster_i} vs {cluster_j}"),
             is_major = transcriptId %in% maj_isoforms$transcript_id ) %>%
      filter(is_major | p_val_adj < 0.05)

    # find genes with >1 transcripts in different clusters
    gene_list <- group_by(markers, geneId) %>%
      summarise(transcripts = n_distinct(transcriptId), clusters = n_distinct(cluster)) %>%
      filter(transcripts >1 & clusters >1)

    # return only markers passing the filters
    output_markers <- semi_join(markers, gene_list, by="geneId")
    return(output_markers)
  }

  cluster_combinations <- as.data.frame(combn(clusters, 2))
  maj_isoforms <- iso_compute_stats(seurat@assays[[assay]]@counts) %>% filter(is_top==TRUE)
  output <- purrr::map_dfr(cluster_combinations, ~ get_markers(.x[1],.x[2]), obj=seurat)
  return(output)
}

# ______________________________________________________________________________

#' Internal helper function to select the criteria used to sort the list of switches
#' centralizing it here, to avoid duplicating code
#'
#' @param switch_list switch list
#' @param order One of "minp_not_major" (default), "max_p", "diff.pct"
#'
#' @return sorted dataframe according to input criteria
#' @export
#' @import tidyr
#'
#' @examples
._order_switches <- function(switch_list, order=NULL){

  if(is.null(order))
    order <- "minp_not_major" # default

  if(order=="max_p")
    return(arrange(switch_list, max_p))
  if (order=="minp_not_major")
    return(arrange(switch_list, minp_not_major))
  if (order=="diff.pct")
    return(arrange(switch_list, desc(diff.pct)))

  stop("check order parameter")
}

# ______________________________________________________________________________


#' Internal worker function for compute_switches
#'
#' @param gene Gene filter
#' @param marker_list input marker list as returned by ISO_SWITCH_ALL
#' @param cutoff unused
#'
#' @return data.frame with isoform switches
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
._gene_switches <- function(gene, marker_list, cutoff=0.05) {

  gene_markers <- marker_list %>% filter(geneId==gene) %>% group_by(contrast) %>% nest()

  output <- data.frame()
  for (i in seq_len(nrow(gene_markers))) {

    row_i <- gene_markers$data[[i]]

    clusters = unique(row_i$cluster)
    group1 <- row_i[row_i$cluster == clusters[1],]$transcriptId
    group2 <- row_i[row_i$cluster == clusters[2],]$transcriptId

    # loop over possible combinations
    for (marker_combination in cross2(group1, group2)) {

      # retrieve full transcript rows
      marker1 <- row_i[row_i$transcriptId == marker_combination[1],]
      marker2 <- row_i[row_i$transcriptId == marker_combination[2],]

      # build switch
      p_val1 = as.numeric(marker1$p_val_adj)
      p_val2 = as.numeric(marker2$p_val_adj)
      log2fc1 = marker1$avg_log2FC
      log2fc2 = marker2$avg_log2FC
      not_major_pval1 = if_else(marker1$is_major, 1.0, p_val1)
      not_major_pval2 = if_else(marker2$is_major, 1.0, p_val2)
      switch <- data.frame(geneId=marker1$geneId,
                           t1=marker1$transcriptId, c1=marker1$cluster, p1=p_val1, log2fc1=log2fc1, m1=marker1$is_major, diff.pct1=marker1$pct.1-marker1$pct.2,
                           t2=marker2$transcriptId, c2=marker2$cluster, p2=p_val2, log2fc2=log2fc2, m2=marker2$is_major, diff.pct2=marker2$pct.1-marker2$pct.2,
                           p_val_adj =  p_val1 * p_val2,
                           max_p = max(p_val1, p_val2),
                           min_log2FC = min(abs(c(log2fc1, log2fc2))),
                           minp_not_major = min(not_major_pval1, not_major_pval2)) %>%
        mutate(diff.pct = abs(diff.pct1-diff.pct2))

      output <- rbind(output , switch)
    }
  }

  if(nrow(output)) {
    output <- ._order_switches(output)
  }

  return(output)
}

# ______________________________________________________________________________

#' Computes isoform switches for a specific gene in a marker list
#'
#' @param marker_list Marker list from ISO_SWITCH_ALL
#' @param cutoff tbc
#' @param gene (optional)
#' @param cluster (optional)
#' @param order (optional)
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
compute_switches <- function(marker_list, cutoff=0.05, gene=NULL, cluster=NULL, order=NULL) {

  if(!is.null(gene))  {
    gene_list<- gene
  } else {
    gene_list <- unique(marker_list$geneId)
  }

  switch_list  <- purrr::map_dfr(gene_list, ._gene_switches, marker_list)

  if(!is.null(cluster))
    switch_list <- filter(switch_list, c1==cluster | c2==cluster)
  if(nrow(switch_list))
    switch_list <- ._order_switches(switch_list, order)

  return(switch_list)
}


# ______________________________________________________________________________

#' Clean up / format switch data frame for reports
#'
#' @param switch_df Output from call to compute_switches
#'
#' @return same df with selected columns and number formatting
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
format_switch_table <- function(switch_df) {
  out <- switch_df %>%
    select(geneId, t1, c1, p1, log2fc1, t2, c2, p2, log2fc2) %>%
    mutate(p1 = format(p1, digits=3),
           log2fc1 = format(log2fc1, digits=2),
           p2 = format(p2, digits=3),
           log2fc2 = format(log2fc2, digits=2))

  return(out)
}




#' Returns a html switch table for reports
#'
#' Can only be used on html output documents.
#'
#' @param marker_list Marker list returned by ISO_SWITCH_ALL
#'
#' @return A reactable printable object
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom reactablefmtr default
#'
#' @examples
gene_switch_table <- function(marker_list) {

  # Find unique genes in marker_list with lowest adj p-value
  react_df <- group_by(marker_list, geneId) %>%
    slice_min(order_by=p_val_adj, n=1, with_ties=FALSE) %>%
    arrange(p_val_adj) %>%
    # filter columns of interest & formatting
    dplyr::select(geneId, avg_log2FC, pct.1, pct.2, p_val_adj, contrast) %>%
    mutate(p_val_adj = as.double(format(p_val_adj, digits=3)))

  table <- htmltools::tagList(
    reactR::html_dependency_react(),
    reactR::html_dependency_reacttools(),
    reactable(react_df,
              style = list(fontFamily = "Work Sans, sans-serif", fontSize = "12px"),
              defaultColDef = colDef(minWidth=80),
              columns = list(
                geneId = colDef(filterable = TRUE, width=90),
                avg_log2FC = colDef(format = colFormat(digits = 2), width=90),
                pct.1 = colDef(width=70),
                pct.2 = colDef(width=70),
                p_val_adj = colDef(width=120, align = "right"),
                contrast = colDef(minWidth = 150, filterable = TRUE)
              ),
              defaultPageSize = 15,
              striped=TRUE,
              compact=TRUE,
              theme = reactablefmtr::default(centered = TRUE, font_size=12, header_font_size = 13))
  )

  return(table)
}

