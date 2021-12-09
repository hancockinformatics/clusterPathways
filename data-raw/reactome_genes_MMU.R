
# Load packages -----------------------------------------------------------

library(tidyverse)


# Load all Reactome pathways/genes ----------------------------------------

reactome_genes_all <- read_tsv(
  "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt",
  col_names = c("gene", "id", "url", "description", "evidence", "species")
)


# Filter for human --------------------------------------------------------

reactome_genes_MMU_0 <- filter(
  reactome_genes_all,
  species == "Mus musculus",
  str_detect(gene, "^ENSMU")
) %>%
  select(-c(url, evidence, species))

reactome_genes_MMU <- reactome_genes_MMU_0 %>%
  group_by(id, description) %>%
  summarise(
    bg_genes = paste(gene, collapse = "; "),
    .groups = "drop"
  )


# Save --------------------------------------------------------------------

usethis::use_data(reactome_genes_MMU, overwrite = TRUE)
