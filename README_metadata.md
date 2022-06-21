
<!-- README.md is generated from README.Rmd. Please edit that file -->

# isoswitch metadata

This page describes how to set set up the metadata required for
isoswitch\_report(), in particular the gtf\_df and transcript\_metadata
parameters

``` r
isoswitch_report(seurat, "isoform", gene="HYAL2", marker_list=switch_markers, gtf_df, transcript_metadata) 
```

## Annotation file (gtf\_df)

isoswitch requires an extract from a standard gene annotation file to
plot gene locus

The GTF raw file can be downloaded from
[GENCODE](https://www.gencodegenes.org/human/)

To optimize file size we’re going to extract the columns we need, and
put them in a dataframe.

``` r
library(rtracklayer)

gtf <- rtracklayer::import('/data/10x_data/10x_5psanger/gencode.v39.annotation.gtf') 
gtf_df <- as.data.frame(gtf) %>%
  dplyr::select(seqnames, start, end, strand, type, gene_name, transcript_id, transcript_name, transcript_type, tag) %>%
  mutate(transcript_id = str_extract(transcript_id, ".*?(?=\\.)"))
```

The data frame can be saved to a local file for convenience:

``` r
local_path <- "/data/10x_data/10x_5psanger/gencode.v39.annotation.rds"
saveRDS(gtf_df, local_path)

# read saved object to use on isoswitch
gtf_df <- readRDS(local_path)
```

## Transcript metadata

isoswitch also needs transcript and gene-level metadata than can be
extracted from ensembl using the `biomart` library

``` r

#  retrieve metadata from biomart to store locally
library(biomaRt)
library(stringr)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# pulls gene names from RNA & isoform assays,
seurat <- readRDS("./seurat_object.rds")
gene_list <- rownames(seurat@assays$RNA)

# Gene-level metadata
gene_metadata <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id', 'description','gene_biotype','transcript_count', 'entrezgene_id'),
                       filters = 'external_gene_name',
                       values = gene_list, 
                       mart = ensembl)

# Transcript metadata
transcript_metadata <- getBM(attributes = c('external_gene_name', 'ensembl_transcript_id','external_transcript_name', 'transcript_biotype','transcript_is_canonical'),
                             filters = 'external_gene_name',
                             values = gene_list,
                             mart = ensembl)

cds_lenghts <- getBM(attributes = c('ensembl_transcript_id','cds_length'),
                             filters = 'external_gene_name',
                             values = gene_list,
                             mart = ensembl)

transcript_metadata <- merge(transcript_metadata, cds_lenghts)

saveRDS(gene_metadata, "./data/gene_metadata.rds")
saveRDS(transcript_metadata, "./data/transcript_metadata.rds")
```
