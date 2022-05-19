#' Compute stats on isoform expression relative to the total expression of a gene.
#' Flags major and top (most expressed) transcripts
#'
#' @param count_matrix Count matrix of seurat object
#'
#' @return data.frame with feature/gene_id/transcript_id, (sum) at transcript level
#'   relative percentage (perc) over total gene expression, boolean flags for major
#'   (maj) and top isoforms
#' @importFrom magrittr %>%
#' @import dplyr
#' @export
#'
#' @examples
iso_compute_stats <- function(count_matrix){

  df <- data.frame(feature = rownames(count_matrix)) %>%
    tidyr::separate(feature, into=c("gene_id", "transcript_id"), sep="\\.\\.", remove=FALSE) %>%
    mutate(sum = Matrix::rowSums(count_matrix))

  if(nrow(filter(df, sum==0)))
    warning("features with 0 expression found in matrix")

  # computes expression at gene-level
  gene_totals <- group_by(df, gene_id) %>%
    dplyr::summarize(total_gene = sum(sum), n_isofs = n(), max_sum = max(sum))

  MAJ_THRESHOLD_MARGIN = 0.4
  df <- left_join(df, gene_totals, by=c("gene_id")) %>%
    mutate(perc = sum / total_gene * 100,
           is_major = perc > (100/n_isofs)*(1-MAJ_THRESHOLD_MARGIN),
           is_top  = sum == max_sum)

  return(df)
}



# ______________________________________________________________________________

#' Retrieve stats (sum, n_cells, avg) on transcripts expressed in a given cell type
#' Transcripts with zero expression are filtered out.
#'
#' @param obj Seurat object
#' @param count_matrix assay count matrix
#' @param cell_type Identity
#' @param features (optional) subset rows in matrix
#'
#' @return data.frame with features expressed in a cell type
#'    (feature/gene_id/transcript_id) + cell_type, sum, n_cells, avg
#' @importFrom magrittr %>%
#' @import dplyr
#' @export
#'
#' @examples
celltype_features <- function(obj, count_matrix, cell_type, features=NULL) {

  cell_list <- Seurat::WhichCells(obj, ident = cell_type)
  cell_type_counts <- if(is.null(features))  count_matrix[, cell_list, drop=FALSE] else count_matrix[features, cell_list, drop=FALSE]
  matches <- stringr::str_match(rownames(cell_type_counts), "(?<gene>.*)\\.\\.(?<transcript>.*)")
  df <- data.frame(row.names=NULL,
                   feature=rownames(cell_type_counts),
                   gene_id=matches[,"gene"],
                   transcript_id=matches[,"transcript"],
                   cell_type=cell_type,
                   sum=Matrix::rowSums(cell_type_counts),
                   n_cells = length(cell_list)) %>%
    mutate(avg = sum/n_cells) %>%
    filter(sum > 0 )
}

# ______________________________________________________________________________

#' Compute number of genes and transcripts per cell type using pseudo-bulk
#'
#' @param obj Seurat object
#' @param assay Assay to pull data from
#' @param bulk_bin_size Size of the pseudobulk groups, defaults to 10.
#'
#' @return Dataframe with number of genes and number of transcripts found in each bin
#' @importFrom magrittr %>%
#' @importFrom Hmisc partition.vector
#' @import dplyr
#' @export
#'
#' @examples
celltype_features_bulk <- function(obj, assay, bulk_bin_size=10) {

  count_matrix <- obj@assays[[assay]]@counts

  # split here all feature names, to avoid redoing it in every iteration
  names_df <- data.frame(feature=rownames(count_matrix)) %>%
    tidyr::separate(feature, into=c("gene_id", "transcript_id"), sep="\\.\\.", remove=FALSE)

  celltype_data = data.frame()

  for(cell_type in levels(obj@active.ident)) {
    # make bins of the same size, dropping remaining (modulo) cells
    cell_list <- Seurat::WhichCells(obj, ident = cell_type)
    n_bins <- floor(length(cell_list) / bulk_bin_size)
    cell_list <- cell_list[1: (n_bins * bulk_bin_size)]
    bins <- partition.vector(cell_list, rep.int(bulk_bin_size, n_bins))

    i = 0
    for(bin in bins) {
      bin_counts <- count_matrix[ ,bin, drop=FALSE]
      expressed_transcripts <- Matrix::rowSums(bin_counts) > 0
      n_genes <- n_distinct(names_df[expressed_transcripts, 'gene_id'])
      n_trans <- n_distinct(names_df[expressed_transcripts, 'transcript_id'])
      # add new row
      new_row <- data.frame(t=cell_type, i=i, n_genes=n_genes, n_transcripts=n_trans)
      celltype_data <- bind_rows(celltype_data, new_row)
      i<-i+1
    }
  }
  return(celltype_data)
}

# ______________________________________________________________________________

#' Clean up / format switch data frame for reports
#'
#' @param switch_df Output from call to compute_switches
#'
#' @return same df with selected columns and number formatting
#' @importFrom magrittr %>%
#' @import dplyr
#' @export
#'
#' @examples
format_switch_table <- function(switch_df) {
  out <- switch_df %>%
    select(t1, c1, p1, log2fc1, t2, c2, p2, log2fc2) %>%
    mutate(p1 = format(p1, digits=3),
           log2fc1 = format(log2fc1, digits=2),
           p2 = format(p2, digits=3),
           log2fc2 = format(log2fc2, digits=2))

  return(out)
}


# gtf_df <- readRDS("/data/10x_data/10x_5psanger/isoswitch/data/gencode.v39.annotation.rds")
# transcript_metadata <- readRDS("/data/10x_data/10x_5psanger/isoswitch/data/adult_iso_metadata.rds")
# test <- readRDS("/data/10x_data/10x_5psanger/output/adult50.rds")
# isoswitch::ISO_SWITCH_ALL(test, clusters=c("Basal", "Suprabasal_1", "EC Venous"), assay="multi")
# isoswitch::plot_marker_matrix(test, mks)
# isoswitch::plot_marker_score(test,mks)
#
