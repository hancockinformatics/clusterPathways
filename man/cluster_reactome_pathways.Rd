% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_reactome_pathways.R
\name{cluster_reactome_pathways}
\alias{cluster_reactome_pathways}
\title{Cluster Reactome Pathways}
\usage{
cluster_reactome_pathways(
  input_pathways,
  input_genes,
  species = "human",
  output_dir = NULL,
  width = 10,
  height = 20
)
}
\arguments{
\item{input_pathways}{Table of Reactome pathways to be clustered. Must have
the pathway ID in the first column, and description in the second. If a
"direction" column (case sensitive) is present, it will be added to heatmap
annotations (see Details for more information).}

\item{input_genes}{Character vector of genes used to generate the
\code{input_pathways} table, e.g. a list of DE genes. Must be Ensembl IDs.}

\item{species}{Either "human" (the default) or "mouse".}

\item{output_dir}{Directory to save heatmaps into. It will be created if it
doesn't already exist.}

\item{width}{Width and height of output heatmaps in inches.}

\item{height}{Width and height of output heatmaps in inches.}
}
\value{
{
A named list containing three data frames, "all_pathways", "rep_pathways",
and "missing_pathways". (see below for more information). The function
also saves heatmaps to ".png" files in the provided output directory.

The additional columns in the first two results tables are:

\describe{
\item{level_1, level_2}{The highest and second-highest levels for each
pathway from the Reactome hierarchy}
\item{n_bg_genes}{The total number of genes annotated to the pathway}
\item{n_cd_genes}{The number genes annotated to the pathway that were
present in the provided input list (candidate genes)}
\item{gene_ratio}{\code{n_cd_genes} / \code{n_bg_genes}}
\item{cluster}{Denotes the group each pathway is placed within}
}

The "rep_pathways" table contains all of the above, plus the column
"n_pathways" which is the number of pathways from that cluster (the number in
square brackets in the "heatmap_rep_pathways_clustered.png" output image).
The "missing_pathways" table contains any input pathways that were not
included in the analysis/figures, typically due to missing data.

}
}
\description{
Using a data frame of Reactome pathways, create a pairwise
Jaccard matrix, then cluster pathways accordingly. Heatmaps/denrograms are
saved, and two results tables are returned as output containing all the
input pathways and which cluster they belonged to.
}
\details{
The direction column must contain either "up" or "down" for each
pathway present. Both column name and contents are case-sensitive.
}
\examples{
\dontrun{
  library(clusterPathways)
  library(tidyverse)

  # Load table of pathways
  my_pathways <- read_csv("pathway_enrichment_result.csv")

  # Load list of genes, which were tested to obtain the above pathways
  my_genes <- read_csv("de_genes.csv") \%>\% pull(ensembl_gene_id)

  # Run cluster_reactome_pathways
  cluster_reactome_pathways(
    input_pathways = my_pathways,
    input_genes = my_genes,
    species = "human",
    output_dir = "clustered_pathways",
    width = 18,
    height = 30
  )
}

}
\references{
None.
}
\seealso{
\url{https://www.github.com/hancockinformatics/clusterPathways}
}
