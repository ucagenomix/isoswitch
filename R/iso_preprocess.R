
# ______________________________________________________________________________

#' Preprocess isoform-level count matrix for isoform switch.
#'
#' Removes transcripts with low relative expression compared to the total gene
#' expression. Filtering threshold is given as parameter. Single-isoform genes
#' (after the filtering step) are also removed. The resulting matrix is stored
#' as a new assay in the input Seurat object
#'
#' @param obj Input Seurat object
#' @param assay Assay containing unfiltered isoform-level counts
#' @param new_assay Name of the new assay to store pre-processed counts
#' @param filter_threshold Filtering threshold expressed as percentage of total gene expression
#'
#' @return Input Seurat object with additional assay contaning pre-processed counts
#' @export
#'
#' @examples
iso_preprocess <- function(obj, assay, new_assay, filter_threshold) {

  count_matrix <- obj@assays[[assay]]@counts
  cat(paste("Input features: ", dim(count_matrix)[[1]]), "\n")
  major_isoforms <- iso_compute_stats(as.data.frame(count_matrix))
  
  # Remove low-expressed features (under threshold param)
  highly_expressed <- filter(major_isoforms, perc > filter_threshold) %>% pull(feature)
  high_expr_matrix <- count_matrix[highly_expressed, , drop=FALSE]

  # Remove features for single-isoform genes
  major_isoforms <- iso_compute_stats(high_expr_matrix)
  multi_features <- filter(major_isoforms, n_isofs > 1) %>% pull(feature)
  only_multis <- high_expr_matrix[multi_features, , drop=FALSE]
  cat(paste("Output dims: ", dim(only_multis)[[1]]), "\n")

  # create new assay with filtered count matrix
  obj[[new_assay]] <- CreateAssayObject(counts = only_multis)

  # Normalize & scale isoform data
  # Use custom scaling for each isoform feature at gene level
  DefaultAssay(obj) <- new_assay
  obj <- NormalizeData(obj)
  obj <- scale_isoform(obj, new_assay)
  return(obj)

}

# ______________________________________________________________________________

#' Scales raw transcript counts relative to total gene expression
#'
#' @param seurat Seurat object
#' @param assay Assay to pull data from
#' @param name Name for new assay with scaled data
#'
#' @return Seurat object with scaled data in new slot
#' @export
#'
#' @author KL
#' @examples
scale_isoform <- function(seurat, assay, name="scale.data"){

  data_in <- as.matrix(seurat@assays[[assay]]@data)

  # retrieve unique gene names as df
  mat <- sapply(strsplit(rownames(data_in), "\\.\\."), `[`, 1)
  all.genes <- unique(mat)
  all.genes <- as.data.frame(all.genes)

  print(paste("Unique features =", dim(all.genes)[1], "out of", dim(data_in)[1], "total features", sep=" "))

  colnames(all.genes)<-c("geneId")
  rownames(all.genes) <- all.genes$geneId

  output <- data_in

  # for all genes
  for(i in 1:dim(all.genes)[1]){

    # print(paste("i",i,sep="="))
    x <- which(mat == all.genes[i,])
    m <- mean(unlist(data_in[x,]))
    sd <- sd(unlist(data_in[x,]))

    # for all gene isoforms
    for(j in 1:length(x)){
      # print(paste("j",j,length(x),sep="="))

      # for all cells
      for(k in 1:dim(data_in)[2]){
        #print(paste("k",k,sep="="))
        output[x[j],k] <- (data_in[x[j],k] - m)/sd
      }
    }
  }

  seurat <- Seurat::SetAssayData(seurat, slot = name, output, assay = assay)
  return(seurat)
}

