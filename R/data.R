#' First two levels of Reactome hierarchy for all pathways in humans
#'
#' A data frame containing all human Reactome pathways, and the corresponding
#' level 1 and 2 names/terms for each.
#'
#' @format A data frame with 2546 rows and 4 columns:
#' \describe{
#'   \item{id}{Reactome ID for all pathways}
#'   \item{description}{Name of the pathway}
#'   \item{level_1}{Highest-level name for the pathway}
#'   \item{level_2}{Second highest-level name for the pathway}
#' }
"reactome_categories_HSA_L1_L2"


#' Genes for all human Reactome pathways
#'
#' A data frame containing the constituent genes for all human Reactome pathways
#'
#' @format A data frame with 2503 rows and 3 columns:
#' \describe{
#'   \item{id}{Reactome ID for all pathways}
#'   \item{description}{Name of the pathway}
#'   \item{bg_genes}{All genes annotated to the pathway, as Ensembl IDs and
#'     separated by "; "}
#' }
"reactome_genes_HSA"


#' First two levels of Reactome hierarchy for all pathways in mice
#'
#' A data frame containing all human Reactome pathways, and the corresponding
#' level 1 and 2 names/terms for each.
#'
#' @format A data frame with 1699 rows and 4 columns:
#' \describe{
#'   \item{id}{Reactome ID for all pathways}
#'   \item{description}{Name of the pathway}
#'   \item{level_1}{Highest-level name for the pathway}
#'   \item{level_2}{Second highest-level name for the pathway}
#' }
"reactome_categories_MMU_L1_L2"


#' Genes for all mouse Reactome pathways
#'
#' A data frame containing the constituent genes for all human Reactome pathways
#'
#' @format A data frame with 1697 rows and 3 columns:
#' \describe{
#'   \item{id}{Reactome ID for all pathways}
#'   \item{description}{Name of the pathway}
#'   \item{bg_genes}{All genes annotated to the pathway, as Ensembl IDs and
#'     separated by "; "}
#' }
"reactome_genes_MMU"
