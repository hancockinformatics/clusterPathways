#' Title
#'
#' @param input_pathways Table of Reactome pathways to be clustered. Must have
#'   the pathway ID in the first column, and description in the second.
#' @param input_genes Character vector of genes used to generate the
#'   `input_pathways` table, e.g. a list of DE genes. Must be Ensembl IDs.
#' @param species Either "human" (the default) or "mouse"
#'
#' @return A data frame (tibble) of gene ratio information for each input
#'   pathway
#'
#' @export
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import purrr
#'
#' @references None.
#'
get_gene_ratio <- function(input_pathways, input_genes, species = "human") {

  if (!is.character(input_genes)) {
    stop("Argument 'input_genes' must be a character vector of Ensembl gene ",
         "IDs.")
  }

  ### Tidy and prep input
  pathway_table_1 <- input_pathways %>%
    remove_rownames() %>%
    rename("id" = 1, "description" = 2)

  if ( any(duplicated(pathway_table_1$id)) ) {
    stop("Your 'input_pathways' contains duplicate IDs. Please remove any ",
         "duplicated pathways and try again.")
  }

  if (species == "human") {
    pathway_table_2 <- pathway_table_1 %>%
      filter(
        id %in% reactome_categories_HSA_L1_L2$id,
        !description %in% reactome_categories_HSA_L1_L2$level_1
      ) %>%
      left_join(reactome_categories_HSA_L1_L2, by = c("id", "description"))


    ### Get the background and candidate genes for each pathway
    pathways_bg_genes <- reactome_genes_HSA %>%
      select(id, bg_genes) %>%
      separate_rows(bg_genes, sep = "; ") %>%
      filter(id %in% pathway_table_2$id) %>%
      split(x = .$bg_genes, f = .$id)

    pathways_cd_genes <- reactome_genes_HSA %>%
      select(id, bg_genes) %>%
      separate_rows(bg_genes, sep = "; ") %>%
      filter(
        id %in% pathway_table_2$id,
        bg_genes %in% input_genes
      ) %>%
      split(x = .$bg_genes, f = .$id)

  } else if (species == "mouse") {
    pathway_table_2 <- pathway_table_1 %>%
      remove_rownames() %>%
      filter(
        id %in% reactome_categories_MMU_L1_L2$id,
        !description %in% reactome_categories_MMU_L1_L2$level_1
      ) %>%
      left_join(reactome_categories_MMU_L1_L2, by = c("id", "description"))


    ### Get the background and candidate genes for each pathway
    pathways_bg_genes <- reactome_genes_MMU %>%
      select(id, bg_genes) %>%
      separate_rows(bg_genes, sep = "; ") %>%
      filter(id %in% pathway_table_2$id) %>%
      split(x = .$bg_genes, f = .$id)

    pathways_cd_genes <- reactome_genes_MMU %>%
      select(id, bg_genes) %>%
      separate_rows(bg_genes, sep = "; ") %>%
      filter(
        id %in% pathway_table_2$id,
        bg_genes %in% input_genes
      ) %>%
      split(x = .$bg_genes, f = .$id)
  } else {
    stop("Argument 'species' must be one of 'human' (default) or 'mouse'.")
  }

  pathway_table_3 <- pathway_table_2 %>%
    mutate(
      n_bg_genes = map_dbl(id, ~length(pathways_bg_genes[[.x]])),
      n_cd_genes = map_dbl(id, ~length(pathways_cd_genes[[.x]])),
      gene_ratio = round((n_cd_genes / n_bg_genes) * 100, digits = 1),
      cd_genes = map_chr(id, ~paste0(pathways_cd_genes[[.x]], collapse = "; ")),
      bg_genes = map_chr(id, ~paste0(pathways_bg_genes[[.x]], collapse = "; "))
    )

  return(pathway_table_3)
}
