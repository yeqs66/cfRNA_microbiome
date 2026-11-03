run_bcaa_network <- function(data, feature_sets, out_dir, repo_root,
                             settings = list()) {
  validate_samples(data)
  allowed <- c("min_n", "max_missing", "min_micro_prevalence", "rho_cutoff",
               "fdr_cutoff", "p_adjust_method", "max_plot_edges", "label_top_n",
               "seed", "include_all")
  if (length(setdiff(names(settings), allowed))) stop("Unknown network settings.")
  defaults <- list(min_n = 20L, max_missing = 0.50, min_micro_prevalence = 0.10,
                   rho_cutoff = 0.10, fdr_cutoff = 0.05, p_adjust_method = "BH",
                   max_plot_edges = 80L, label_top_n = 20L, seed = 2026L,
                   include_all = TRUE)
  defaults[names(settings)] <- settings
  env <- source_functions(file.path(repo_root, "reference/bcaa_network/02.R"),
                          "run_cross_omics_network")
  do.call(env$run_cross_omics_network,
          c(list(matrix_meta = data, feature_sets = feature_sets, out_dir = out_dir), defaults))
}
