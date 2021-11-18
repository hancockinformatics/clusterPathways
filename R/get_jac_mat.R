#' Create matrix for Jaccard calculation
#'
#' @param list A named list of pathways and their constituent genes
#'
#' @import purrr
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @description An internal helper function to construct a matrix from the list
#' of pathways and their genes.
#'
#' @return A matrix of overlapping genes among all pathways
#'
#' @references None.
#'
#' @seealso <https://www.github.com/hancockinformatics/clusterReactome>
#'
get_jac_mat <- function(list) {
  list %>%
    map(~data.frame(id = .x)) %>%
    bind_rows(.id = "name") %>%
    mutate(present = 1) %>%
    pivot_wider(
      id_cols     = "id",
      names_from  = "name",
      values_from = "present"
    ) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "id")
}
