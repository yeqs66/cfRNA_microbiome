# Infrastructure only. Statistical formulas are kept in source-backed modules.
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

validate_samples <- function(d, key = "sign") {
  if (!is.data.frame(d) || !key %in% names(d) || anyDuplicated(names(d)))
    stop("Expected a data.frame with unique columns and key: ", key)
  k <- as.character(d[[key]])
  if (anyNA(k) || any(!nzchar(k)) || anyDuplicated(k))
    stop("Missing/duplicated sample keys; repair upstream joins, do not drop duplicates arbitrarily.")
  if (key == "sign" && any(!grepl("^.+_(early|middle|late)$", k)))
    stop("sign must end in _early, _middle or _late.")
  invisible(d)
}

read_frame <- function(path) {
  if (is.null(path) || length(path) != 1L || !file.exists(path))
    stop("Input unavailable: ", paste(path, collapse = ", "))
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

source_functions <- function(path, symbols, env = new.env(parent = globalenv())) {
  expressions <- parse(path, encoding = "UTF-8")
  found <- character()
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% symbols &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
      eval(expr, env)
      found <- c(found, as.character(expr[[2]]))
    }
  }
  if (length(setdiff(symbols, found))) stop("Source functions missing: ",
    paste(setdiff(symbols, found), collapse = ", "))
  env
}

write_run_metadata <- function(out_dir, config, stage, mode, status) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(stage = stage, mode = mode, status = status, config = config,
               time = format(Sys.time(), tz = "UTC", usetz = TRUE)),
          file.path(out_dir, "run_metadata.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
}

resolve_inputs <- function(x, base) {
  if (is.list(x)) return(lapply(x, resolve_inputs, base = base))
  if (is.character(x)) return(vapply(x, function(p) {
    if (grepl("^(?:[A-Za-z]:[/\\\\]|/|~)", p, perl = TRUE)) p else file.path(base, p)
  }, character(1), USE.NAMES = FALSE))
  x
}
