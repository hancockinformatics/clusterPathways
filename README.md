# clusterPathways
Use Jaccard index and constituent genes to identify groups of similar pathways

## Installation
```r
remotes::install_github("hancockinformatics/clusterPathways")
```

## Example
```r
library(clusterPathways)
library(tidyverse)


# Load table of pathways
my_pathways <- read_csv("pathway_enrichment_result.csv")

glimpse(my_pathways)
#> Rows: 192
#> Columns: 5
#> $ id          <chr> "R-HSA-156842", "R-HSA-156902", "R-HSA-192823", "R-HSA-9…
#> $ description <chr> "Eukaryotic Translation Elongation", "Peptide chain elon…
#> $ direction   <chr> "up", "up", "up", "up", "up", "up", "up", "up", "up", "u…
#> $ pvalue      <dbl> 1.400442e-39, 4.774000e-38, 6.119370e-35, 7.946710e-35, …
#> $ p_adjust    <dbl> 1.819174e-36, 3.100713e-35, 2.580694e-32, 2.580694e-32, …


# Load list of genes, which were tested to obtain the above pathways
my_genes <- read_csv("de_genes.csv") %>% pull(ensembl_gene_id)
head(my_genes)
#> [1] "ENSG00000100664" "ENSG00000114529" "ENSG00000055147" "ENSG00000186226" 
#> [5] "ENSG00000243349" "ENSG00000144136"


# Run cluster_reactome_pathways
cluster_reactome_pathways(
  input_pathways = my_pathways,
  input_genes    = my_genes,
  species        = "human",
  output_dir     = "clustered_pathways",
  width  = 18,
  height = 30
)
```
