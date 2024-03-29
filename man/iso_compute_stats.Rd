% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iso_stats.R
\name{iso_compute_stats}
\alias{iso_compute_stats}
\title{Compute stats on isoform expression relative to the total expression of a gene.
Flags major and top (most expressed) transcripts}
\usage{
iso_compute_stats(count_matrix, count_zeros = TRUE)
}
\arguments{
\item{count_matrix}{Count matrix of seurat object}

\item{count_zeros}{if TRUE (default), features with zero expression are included in the stats}
}
\value{
data.frame with feature/gene_id/transcript_id, (sum) at transcript level
relative percentage (perc) over total gene expression, boolean flags for major
(maj) and top isoforms
}
\description{
Compute stats on isoform expression relative to the total expression of a gene.
Flags major and top (most expressed) transcripts
}
