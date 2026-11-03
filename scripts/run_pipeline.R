# Usage: Rscript scripts/run_pipeline.R --config config/demo.yml --stage network --mode demo
if (.Platform$OS.type == "windows") {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE))
}
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("No analysis runs by default. Specify --config PATH --stage NAME --mode demo|real explicitly.")
get_arg <- function(key, default = NULL) {
  idx <- which(args == key)
  if (!length(idx)) return(default)
  if (length(idx) != 1L || idx == length(args)) stop("Invalid argument: ", key)
  args[idx + 1L]
}
unknown <- args[seq.int(1L, length(args), by = 2L)]
if (length(args) %% 2L || any(!unknown %in% c("--config", "--stage", "--mode")))
  stop("Arguments must be --config PATH --stage NAME --mode demo|real")
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."),
                           winslash = "/", mustWork = TRUE)
for (f in c("core.R", "network.R", "metabolomics.R", "limma.R", "demo.R", "reference_adapter.R"))
  source(file.path(repo_root, "R", f), encoding = "UTF-8")
require_packages("yaml")
config_path <- normalizePath(get_arg("--config", file.path(repo_root, "config/demo.yml")),
                             winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(config_path)
stage <- get_arg("--stage")
mode <- get_arg("--mode")
if (is.null(stage) || is.null(mode)) stop("Explicit --stage and --mode required.")
if (!mode %in% c("demo", "real")) stop("mode must be demo or real")
if (is.null(config$output_dir)) stop("output_dir required")
config_base <- dirname(config_path)
output_root <- resolve_inputs(config$output_dir, config_base)
out_dir <- file.path(output_root, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"),
                                        "_", Sys.getpid(), "_", mode, "_", stage))
if (!grepl("^[A-Za-z0-9_-]+$", stage)) stop("Invalid stage name")
if (dir.exists(out_dir)) stop("Run directory already exists; wait or choose another output directory")
demo <- if (mode == "demo") make_demo() else NULL
if (mode == "real" && !is.null(config$inputs)) config$inputs <- resolve_inputs(config$inputs, config_base)
status <- "failed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tryCatch({
  if (stage == "network") {
    data <- if (mode == "demo") demo$network else read_frame(config$inputs$network)
    fs <- if (mode == "demo") demo$feature_sets else config$feature_sets
    result <- run_bcaa_network(data, fs, out_dir, repo_root,
                               if (is.null(config$network)) list() else config$network)
  } else if (stage == "metabolomics_lmm") {
    data <- if (mode == "demo") demo$lmm else read_frame(config$inputs$metabolomics)
    map <- if (mode == "demo") demo$metabolite_map else read_frame(config$inputs$metabolite_map)
    result <- run_metabolite_lmm(data, map, out_dir, repo_root)
  } else if (stage %in% c("limma_mirna", "limma_lncrna")) {
    if (mode == "real" && !isTRUE(config$acknowledge_legacy_limma))
      stop("Read audit L02-L04 first; historical contrasts/filter are preserved. Explicit acknowledge_legacy_limma: true required.")
    counts <- if (mode == "demo") demo$counts else readRDS(config$inputs$counts)
    metadata <- if (mode == "demo") demo$metadata else read_frame(config$inputs$metadata)
    result <- run_limma_reference(counts, metadata, sub("limma_", "", stage), out_dir)
  } else if (stage %in% names(config$reference_stages)) {
    if (mode == "demo") stop("This reference stage has no validated simulated fixture; see coverage table.")
    recipe <- config$reference_stages[[stage]]
    if (!is.null(recipe$paths)) recipe$paths <- resolve_inputs(recipe$paths, config_base)
    if (!is.null(recipe$globals_rds)) recipe$globals_rds <- resolve_inputs(recipe$globals_rds, config_base)
    result <- run_reviewed_reference(recipe, repo_root, out_dir)
  } else {
    stop("Stage not implemented/confirmed: ", stage,
      ". Available standalone: network, metabolomics_lmm, limma_mirna, limma_lncrna. Other stages require an explicitly reviewed reference recipe; see audit.")
  }
  status <- "completed"
}, finally = {
  write_run_metadata(out_dir, config, stage, mode, status)
  message("Run status: ", status, "; output: ", out_dir)
})
