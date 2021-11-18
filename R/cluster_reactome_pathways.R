#' Cluster Reactome pathways using the Jaccard index
#'
#' @param input_df Unaltered data frame of ReactomePA results. If there is a
#'   "direction" column with up/down for each pathway, this information will be
#'   added to the heatmaps.
#' @param reactome_levels Data frame containing the first two levels for all
#'   Reactome pathways in humans. Must have columns "id", "description",
#'   "level_1", "level_2". If you don't have this information, the description
#'   section contains a link to a script that will create this table for you.
#' @param output_dir Directory to save results. It will be created if it doesn't
#'   already exist.
#' @param width Width of main output file in inches.
#' @param height Height of main output file in inches.
#'
#' @export
#'
#' @import pheatmap
#' @import tibble
#' @import dplyr
#' @import stringr
#' @import purrr
#' @import vegan
#' @import janitor
#' @import glue
#' @import RColorBrewer
#'
#' @description Using the data frame of results output from
#'   [ReactomePA](https://www.bioconductor.org/packages/ReactomePA/), create a
#'   pairwise Jaccard matrix and cluster pathways accordingly.
#'   Heatmaps/denrograms are saved, and two results tables are returned as
#'   output.
#'
#' @return {
#'  A named list containing two data frames, "all_pathways" and "rep_pathways"
#'  (see Details for more information). The function also saves heatmaps to
#'  ".png" files in the provided output directory.
#'
#' \describe{
#'   \item{cluster}{Denotes the groups the pathways are placed within}
#'   \item{level_1, level_2}{The highest and second-highest levels for each
#'     pathway from the Reactoem hierarchy}
#'   \item{candidate_genes}{The number genes annotated to the pathway that were
#'     present in the list provided as input to ReactomePA}
#'   \item{genes_in_pathway}{The total number of genes annotated to the pathway}
#'   \item{gene_ratio}{`candidate_genes` / `genes_in_pathway`}
#' }
#'
#' The "rep_pathways" table contains all of the above, plus the column
#' "n_pathways" which is the number of pathways from that cluster (the number in
#' square brackets in the "heatmap_rep_pathways_clustered.png" output image).
#'
#' }
#'
#' @details The following columns must be present in "input_df": "ID",
#'   "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue",
#'   "geneID", "Count".
#'
#'   You will likely need to play with different values of width/height to get
#'   images that are the proper size. In some cases it may be easier to modify
#'   the png() calls directly, and ignore the arguments here.
#'
#'   Link to the script to create the Reactome levels table:
#'   <https://github.com/hancockinformatics/misc_R_scripts/blob/master/R_scripts/get_reactome_L1_L2.R>
#'
#' @references None.
#'
#' @seealso <https://www.github.com/hancockinformatics/clusterReactome>
#'
cluster_reactome_pathways <- function(
  input_df,
  reactome_levels,
  output_dir,
  width,
  height
) {

  # Check inputs ----------------------------------------------------------

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  colnames_to_check <- c(
    "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue",
    "geneID", "Count"
  )

  if (!all(colnames_to_check %in% colnames(input_df))) {
    stop(
      "Argument 'input_df' must contain the following columns: ",
      paste0(colnames_to_check, collapse = ", "),
      "."
    )
  }

  if (!all(
    c("id", "description", "level_1", "level_2") %in% colnames(reactome_levels)
  )) {
    stop(paste0(
      "Argument 'reactome_levels' must contain the following columns: ",
      "'id', 'level_1,' and 'level_2.'"
    ))
  }

  message("Inputs OK...\n")


  # Define helper functions -----------------------------------------------

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


  # Clean inputs ----------------------------------------------------------

  message("Cleaning inputs...")
  input_df_clean <- input_df %>%
    remove_rownames() %>%
    clean_names() %>%
    as.data.frame() %>%
    filter(!description %in% reactome_levels$level_1) %>%
    left_join(select(reactome_levels, -description), by = "id")

  input_list <- input_df_clean %>%
    select(id, gene_id) %>%
    deframe() %>%
    str_split(., "/") %>%
    set_names(input_df_clean$id)

  n_level2 <- input_df_clean %>%
    distinct(level_2) %>%
    nrow()

  message("Done.\n")


  # Create jaccard matrix -------------------------------------------------

  jaccard_mat_id <- get_jac_mat(input_list) %>%
    t() %>%
    vegdist(method = "jaccard", binary = TRUE, diag = TRUE) %>%
    as.matrix()


  # Convert IDs in matrix to pathway names --------------------------------

  jaccard_mat_description <- jaccard_mat_id %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "id_2", values_to = "dist") %>%
    left_join(select(input_df_clean, id, description), by = "id") %>%
    left_join(select(input_df_clean, "id_2" = id, description), by = "id_2") %>%
    pivot_wider(
      description.x,
      names_from = description.y,
      values_from = "dist"
    ) %>%
    column_to_rownames("description.x")


  ### Check for a direction column, and set up annotations if present
  if ("direction" %in% colnames(input_df)) {
    message("Found direction column to be added to the heatmaps...\n")
    ann_colour_table <- input_df_clean %>%
      select(description, direction) %>%
      column_to_rownames("description")
    ann_colour_list <-
      list(direction = c("down" = "springgreen3", "up" = "firebrick"))
  } else {
    message("No direction column identified...\n")
    ann_colour_table <- NULL
    ann_colour_list  <- NULL
  }


  # Create the initial heatmap --------------------------------------------

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


  # Find representative pathways ------------------------------------------

  message("Finding representative pathways and plotting...")
  initial_clusters <- heatmap_initial$tree_row %>%
    cutree(k = n_level2) %>%
    as.data.frame() %>%
    set_names("cluster") %>%
    rownames_to_column("description")

  pathways_w_GR <- input_df_clean %>%
    select(
      id,
      description,
      "candidate_genes" = count,
      bg_ratio,
      level_1,
      level_2,
      pvalue,
      p_adjust
    ) %>%
    mutate(
      genes_in_pathway = as.numeric(str_remove(bg_ratio, "/[0-9]{1,5}$")),
      gene_ratio = candidate_genes / genes_in_pathway
    ) %>%
    select(-bg_ratio)

  initial_clusters_max_GR <- pathways_w_GR %>%
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

  rep_pathway_mat <- jaccard_mat_description[
    initial_clusters_max_GR_chr$description,
    initial_clusters_max_GR_chr$description
  ]

  rep_pathway_mat_n <- rep_pathway_mat
  colnames(rep_pathway_mat_n) <- initial_clusters_max_GR_chr$description_n
  rownames(rep_pathway_mat_n) <- initial_clusters_max_GR_chr$description_n

  if (!is.null(ann_colour_table)) {
    ann_colour_table_n <- input_df_clean %>%
      left_join(initial_clusters_max_GR_chr, by = "description") %>%
      select(description_n, direction) %>%
      filter(!is.na(description_n)) %>%
      column_to_rownames("description_n")
  }

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
    height = height * 0.6,
    units  = "in",
    res    = 150
  )
  print(heatmap_reps)
  dev.off()
  message("Done.\n")


  # Create a heatmap per level_1 group ------------------------------------

  message("Creating heatmaps per Level 1 terms...")
  level_1_groups <- input_df_clean %>%
    select(description, level_1) %>%
    split(x = .$description, f = .$level_1) %>%
    keep(~length(.x) >= 5)

  level_1_mats <- level_1_groups %>% map(function(x) {
    jaccard_mat_description[x, x]
  })

  level_1_n_clust <- input_df_clean %>%
    select(level_1, level_2) %>%
    split(x = .$level_2, f = .$level_1) %>%
    map(~length(unique(.x))) %>%
    magrittr::extract(names(level_1_groups))

  pwalk(
    list(level_1_mats, names(level_1_mats), level_1_n_clust),
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
    }
  )

  ### Leftover pathways in one group
  leftover_pathways <- input_df_clean %>%
    select(description, level_1) %>%
    split(x = .$description, f = .$level_1) %>%
    keep(~length(.x) < 5) %>%
    flatten() %>%
    as.character()

  leftover_sig_jac_mat <-
    jaccard_mat_description[leftover_pathways, leftover_pathways]

  leftover_n_clust <- input_df_clean %>%
    filter(description %in% leftover_pathways) %>%
    pull(level_2) %>%
    unique() %>%
    length()

  png(
    glue("{output_dir}/level1_heatmap_miscellaneous.png"),
    width  = 16,
    height = 12,
    units  = "in",
    res    = 150
  )
  pheatmap(
    mat                  = leftover_sig_jac_mat,
    show_colnames        = FALSE,
    color                = heatmaps_colours,
    annotation_row       = ann_colour_table,
    annotation_colors    = ann_colour_list,
    annotation_names_row = FALSE,
    fontsize             = 14,
    treeheight_col       = 0,
    cutree_rows          = leftover_n_clust,
    main                 = "Miscellaneous pathways"
  )
  dev.off()
  message("Done.\n")


  # Create output tables --------------------------------------------------

  message("Generate output tables...")
  table_1 <- pathways_w_GR %>%
    left_join(initial_clusters, by = "description") %>%
    arrange(cluster, desc(gene_ratio)) %>%
    relocate(
      cluster,
      description,
      id,
      level_1,
      level_2,
      candidate_genes,
      genes_in_pathway,
      gene_ratio,
      pvalue,
      p_adjust
    )

  table_2 <- initial_clusters_max_GR %>%
    select(-description_n) %>%
    filter(description %in% initial_clusters_max_GR_chr$description) %>%
    relocate(
      cluster,
      n_pathways,
      description,
      id,
      level_1,
      level_2,
      candidate_genes,
      genes_in_pathway,
      gene_ratio,
      pvalue,
      p_adjust
    )
  message("Done.")

  return(list(
    all_pathways = table_1,
    rep_pathways = table_2
  ))
}
