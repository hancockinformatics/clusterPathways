.onAttach <- function(...) {

  packageStartupMessage(paste0(
    "Thanks for using clusterPathways v",
    utils::packageVersion("clusterPathways"), "!\n",
    "If you encounter any bugs or problems, please submit an issue at the\n",
    "Github page: https://github.com/hancockinformatics/clusterPathways"
  ))
}
