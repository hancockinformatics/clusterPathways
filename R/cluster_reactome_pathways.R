
#' Cluster Reactome Pathways
#'
#' @param input_pathways Table of Reactome pathways to be clustered. Must have
#'   the pathway ID in the first column, and description in the second. If a
#'   "direction" column (case sensitive) is present, it will be added to heatmap
#'   annotations (see Details for more information).
#' @param input_genes Character vector of genes used to generate the
#'   `input_pathways` table, e.g. a list of DE genes. Must be Ensembl IDs.
#' @param species Either "human" (the default) or "mouse".
#' @param output_dir Directory to save heatmaps into. It will be created if it
#'   doesn't already exist.
#' @param width Width and height of output heatmaps in inches.
#' @param height Width and height of output heatmaps in inches.
#'
#' @description Using a data frame of Reactome pathways, create a pairwise
#'   Jaccard matrix, then cluster pathways accordingly. Heatmaps/denrograms are
#'   saved, and two results tables are returned as output containing all the
#'   input pathways and which cluster they belonged to.
#'
#' @export
#'
#' @import dplyr
#' @import glue
#' @import janitor
#' @import pheatmap
#' @import purrr
#' @import RColorBrewer
#' @import rmarkdown
#' @import stringr
#' @import tibble
#' @import vegan
#'
#' @return {
#'  A named list containing three data frames, "all_pathways", "rep_pathways",
#'  and "missing_pathways". (see below for more information). The function
#'  also saves heatmaps to ".png" files in the provided output directory.
#'
#'  The additional columns in the first two results tables are:
#'
#' \describe{
#'   \item{level_1, level_2}{The highest and second-highest levels for each
#'     pathway from the Reactome hierarchy}
#'   \item{n_bg_genes}{The total number of genes annotated to the pathway}
#'   \item{n_cd_genes}{The number genes annotated to the pathway that were
#'     present in the provided input list (candidate genes)}
#'   \item{gene_ratio}{`n_cd_genes` / `n_bg_genes`}
#'   \item{cluster}{Denotes the group each pathway is placed within}
#' }
#'
#' The "rep_pathways" table contains all of the above, plus the column
#' "n_pathways" which is the number of pathways from that cluster (the number in
#' square brackets in the "heatmap_rep_pathways_clustered.png" output image).
#' The "missing_pathways" table contains any input pathways that were not
#' included in the analysis/figures, typically due to missing data.
#'
#' }
#'
#' @details The direction column must contain either "up" or "down" for each
#'   pathway present. Both column name and contents are case-sensitive.
#'
#' @references None.
#'
#' @examples
#' \dontrun{
#'   library(clusterPathways)
#'   library(tidyverse)
#'
#'   # Load table of pathways
#'   my_pathways <- read_csv("pathway_enrichment_result.csv")
#'
#'   # Load list of genes, which were tested to obtain the above pathways
#'   my_genes <- read_csv("de_genes.csv") %>% pull(ensembl_gene_id)
#'
#'   # Run cluster_reactome_pathways
#'   cluster_reactome_pathways(
#'     input_pathways = my_pathways,
#'     input_genes = my_genes,
#'     species = "human",
#'     output_dir = "clustered_pathways",
#'     width = 18,
#'     height = 30
#'   )
#' }
#'
#' @seealso <https://www.github.com/hancockinformatics/clusterPathways>
#'
cluster_reactome_pathways <- function(input_pathways,
                                      input_genes,
                                      species = "human",
                                      output_dir = NULL,
                                      width = 10,
                                      height = 20) {


  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
  } else {
    stop("Must specify an output directory for results.")
  }

  if (!is.character(input_genes)) {
    stop("Argument 'input_genes' must be a character vector of Ensembl gene ",
         "IDs.")
  }


  ### Tidy and prep input
  message("Tidying input...")

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
      rename("id" = 1, "description" = 2) %>%
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


  ### Find number of level 2 terms represented in input pathways
  n_level2 <- length(unique(pathway_table_2$level_2))


  ### Combine above info into a single table
  pathway_table_3 <- pathway_table_2 %>%
    mutate(
      n_bg_genes = map_dbl(id, ~length(pathways_bg_genes[[.x]])),
      n_cd_genes = map_dbl(id, ~length(pathways_cd_genes[[.x]])),
      gene_ratio = round((n_cd_genes / n_bg_genes) * 100, digits = 1)
    )
  message("Done.\n")


  ### Create Jaccard matrix using internal function and `vegan` function
  jaccard_mat_id <- pathways_cd_genes %>%
    get_jac_mat() %>%
    t() %>%
    vegan::vegdist(
      method = "jaccard",
      binary = TRUE,
      diag   = TRUE
    ) %>%
    as.matrix()


  ### Convert the row/column names of matrix from pathway ID to description
  jaccard_mat_description <- jaccard_mat_id %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "id_2", values_to = "dist") %>%
    left_join(select(pathway_table_2, id, description), by = "id") %>%
    left_join(select(pathway_table_2, "id_2" = id, description), by = "id_2") %>%
    pivot_wider(
      description.x,
      names_from = description.y,
      values_from = "dist"
    ) %>%
    column_to_rownames("description.x")


  ### If present, set up the direction heatmap annotation
  if ("direction" %in% colnames(pathway_table_2)) {
    message("Found direction column to be added to the heatmaps...\n")
    ann_colour_table <- pathway_table_2 %>%
      select(description, direction) %>%
      mutate(direction = str_to_lower(direction)) %>%
      column_to_rownames("description")
    ann_colour_list <-
      list(direction = c("down" = "springgreen3", "up" = "firebrick"))
  } else {
    message("No direction column identified...\n")
    ann_colour_table <- NULL
    ann_colour_list  <- NULL
  }


  ### Create the initial heatmap of all pathways
  message("Creating first heatmap...")
  heatmaps_colours <-
    colorRampPalette(brewer.pal(n = 9, name = "Blues"))(10)

  heatmap_initial <- pheatmap(
    mat                  = jaccard_mat_description,
    show_colnames        = FALSE,
    color                = heatmaps_colours,
    annotation_row       = ann_colour_table,
    annotation_colors    = ann_colour_list,
    annotation_names_row = FALSE,
    fontsize             = 14,
    cutree_rows          = n_level2,
    treeheight_col       = 0,
    silent               = TRUE
  )

  png(
    glue("{output_dir}/heatmap_all_pathways_clustered.png"),
    width  = width,
    height = height,
    units  = "in",
    res    = 150
  )
  print(heatmap_initial)
  dev.off()
  message("Done.\n")


  ### Retrieve clusters from initial heatmap
  message("Finding representative pathways and plotting...")
  initial_clusters <- heatmap_initial$tree_row %>%
    cutree(k = n_level2) %>%
    as.data.frame() %>%
    set_names("cluster") %>%
    rownames_to_column("description")


  ### Find the pathway for each cluster with the highest gene ratio
  initial_clusters_max_GR <- pathway_table_3 %>%
    left_join(initial_clusters, by = "description") %>%
    group_by(cluster) %>%
    mutate(
      n_pathways = n(),
      description_n = case_when(
        n_pathways != 1 ~ glue("{description} [{n_pathways}]"),
        TRUE ~ description
      )
    ) %>%
    ungroup()

  initial_clusters_max_GR_chr <- initial_clusters_max_GR %>%
    group_by(cluster) %>%
    arrange(desc(gene_ratio)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(description, description_n)


  ### Subset the matrix to just the representative pathways identified above
  rep_pathway_mat <- jaccard_mat_description[
    initial_clusters_max_GR_chr$description,
    initial_clusters_max_GR_chr$description
  ]

  rep_pathway_mat_n <- rep_pathway_mat
  colnames(rep_pathway_mat_n) <-
    initial_clusters_max_GR_chr$description_n
  rownames(rep_pathway_mat_n) <-
    initial_clusters_max_GR_chr$description_n

  if (!is.null(ann_colour_table)) {
    ann_colour_table_n <- pathway_table_2 %>%
      left_join(initial_clusters_max_GR_chr, by = "description") %>%
      select(description_n, direction) %>%
      filter(!is.na(description_n)) %>%
      column_to_rownames("description_n")
  } else {
    ann_colour_table_n <- NULL
  }


  ### Create the representative heatmap
  heatmap_reps <- pheatmap(
    mat                  = rep_pathway_mat_n,
    show_colnames        = FALSE,
    color                = heatmaps_colours,
    annotation_row       = ann_colour_table_n,
    annotation_colors    = ann_colour_list,
    annotation_names_row = FALSE,
    fontsize             = 14,
    cutree_rows          = n_level2,
    treeheight_col       = 0,
    silent               = TRUE
  )

  png(
    glue("{output_dir}/heatmap_rep_pathways_clustered.png"),
    width  = width,
    height = height,
    units  = "in",
    res    = 150
  )
  print(heatmap_reps)
  dev.off()
  message("Done.\n")


  ### Create heatmaps per Level 1 pathway/group
  message("Creating heatmaps per Level 1 terms...")
  level_1_groups <- pathway_table_2 %>%
    select(description, level_1) %>%
    split(x = .$description, f = .$level_1) %>%
    keep( ~ length(.x) >= 5)

  level_1_mats <- level_1_groups %>% map(function(x) {
    jaccard_mat_description[x, x]
  })

  level_1_n_clust <- pathway_table_2 %>%
    select(level_1, level_2) %>%
    split(x = .$level_2, f = .$level_1) %>%
    map( ~ length(unique(.x))) %>%
    magrittr::extract(names(level_1_groups))

  pwalk(list(level_1_mats, names(level_1_mats), level_1_n_clust),
        function(input_mat, mat_name, mat_n_clust) {
          png(
            glue("{output_dir}/level1_heatmap_{str_replace_all(mat_name, ' ', '_')}.png"),
            width  = 16,
            height = 12,
            units  = "in",
            res    = 150
          )
          pheatmap(
            mat                  = input_mat,
            show_colnames        = FALSE,
            color                = heatmaps_colours,
            annotation_row       = ann_colour_table,
            annotation_colors    = ann_colour_list,
            annotation_names_row = FALSE,
            fontsize             = 14,
            treeheight_col       = 0,
            cutree_rows          = mat_n_clust,
            main                 = mat_name
          )
          dev.off()
        })

  ### Put all the "leftover" pathways in one group/heatmap
  leftover_pathways <- pathway_table_2 %>%
    select(description, level_1) %>%
    split(x = .$description, f = .$level_1) %>%
    keep( ~ length(.x) < 5) %>%
    flatten() %>%
    as.character()

  leftover_sig_jac_mat <-
    jaccard_mat_description[leftover_pathways, leftover_pathways]

  leftover_n_clust <- pathway_table_2 %>%
    filter(description %in% leftover_pathways) %>%
    pull(level_2) %>%
    unique() %>%
    length()

  heatmap_leftover <- pheatmap(
    mat                  = leftover_sig_jac_mat,
    show_colnames        = FALSE,
    color                = heatmaps_colours,
    annotation_row       = ann_colour_table,
    annotation_colors    = ann_colour_list,
    annotation_names_row = FALSE,
    fontsize             = 14,
    treeheight_col       = 0,
    cutree_rows          = leftover_n_clust,
    main                 = "Miscellaneous pathways",
    silent               = TRUE
  )
  png(
    glue("{output_dir}/level1_heatmap_miscellaneous.png"),
    width  = 16,
    height = 12,
    units  = "in",
    res    = 150
  )
  print(heatmap_leftover)
  dev.off()
  message("Done.\n")


  ## Create output tables and summary Rmd report
  message("Generating output tables...")
  out_table_1 <- initial_clusters_max_GR %>%
    select(-description_n) %>%
    mutate(
      cd_genes = map(id, ~paste0(pathways_cd_genes[[.x]], collapse = "; ")),
      bg_genes = map(id, ~paste0(pathways_bg_genes[[.x]], collapse = "; "))
    )

  out_table_2 <- out_table_1 %>%
    filter(description %in% initial_clusters_max_GR_chr$description)

  out_table_3 <- pathway_table_1 %>%
    filter(!id %in% out_table_1$id)

  rmarkdown::render(
    input = system.file("rmd", "results.Rmd", package = "clusterPathways"),
    params = list(
      table_1 = select(out_table_1, -c(cd_genes, bg_genes)),
      table_2 = select(out_table_2, -c(cd_genes, bg_genes))
    ),
    output_dir = output_dir
  )

  message("Done.\n")
  return(list(
    all_pathways = out_table_1,
    rep_pathways = out_table_2,
    missing_pathways = out_table_3
  ))
}
