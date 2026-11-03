# Explicit, opt-in adapter for historical scripts after the user reviews the audit.
# Does not infer missing models or change statistical expressions.
run_reviewed_reference <- function(recipe, repo_root, out_dir) {
  if (!isTRUE(recipe$reviewed)) stop("Reference stage awaits manual audit. Set reviewed: true only after checking docs/manuscript_v9_code_audit.md and the input recipe.")
  path <- normalizePath(file.path(repo_root, recipe$reference), winslash = "/", mustWork = TRUE)
  reference_root <- paste0(normalizePath(file.path(repo_root, "reference"), winslash = "/"), "/")
  if (!startsWith(path, reference_root)) stop("Reference must stay inside reference/.")
  ex <- parse(path, encoding = "UTF-8")
  mapping <- recipe$paths
  if (is.null(mapping)) mapping <- list()
  env <- new.env(parent = globalenv())
  if (!is.null(recipe$globals_rds)) {
    objects <- readRDS(recipe$globals_rds)
    if (!is.list(objects) || is.null(names(objects)) || any(!nzchar(names(objects))) || anyDuplicated(names(objects)))
      stop("globals_rds must contain a uniquely named list of required upstream objects.")
    list2env(objects, env)
  }
  rewrite <- function(node) {
    if (is.character(node)) {
      for (i in seq_along(node)) {
        if (node[i] %in% names(mapping)) node[i] <- mapping[[node[i]]]
        else if (grepl("__PATH_[0-9]+__|__LOCAL_PATH_REDACTED__", node[i]))
          stop("Unconfigured reference path token: ", node[i])
      }
      return(node)
    }
    if (is.call(node)) {
      head <- as.character(node[[1]])
      if (length(head) == 1L && head == "setwd") return(quote(invisible(NULL)))
      if (length(head) == 1L && head == "rm" && "list" %in% names(as.list(node)))
        return(quote(invisible(NULL)))
      # save.image() captures global session, not the isolated script environment.
      # Suppress this incidental workspace dump; statistical result exports remain.
      if (length(head) == 1L && head == "save.image") return(quote(invisible(NULL)))
      for (i in seq_along(node)) {
        if (identical(node[[i]], quote(expr = ))) next
        node[[i]] <- rewrite(node[[i]])
      }
      return(node)
    }
    if (is.expression(node) || is.pairlist(node)) {
      for (i in seq_along(node)) {
        if (identical(node[[i]], quote(expr = ))) next
        node[[i]] <- rewrite(node[[i]])
      }
    }
    node
  }
  ex <- rewrite(ex) # preflight every path token BEFORE executing any expression
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old <- getwd(); on.exit(setwd(old), add = TRUE)
  setwd(out_dir)
  for (expr in ex) eval(expr, env)
  invisible(TRUE)
}
