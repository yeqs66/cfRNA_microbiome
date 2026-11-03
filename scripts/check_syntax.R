# Static check ONLY: parse code, never source/evaluate scripts or load study data.
if (.Platform$OS.type == "windows") {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE))
}
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), winslash = "/")
files <- c(list.files(file.path(root, "R"), pattern = "\\.[Rr]$", full.names = TRUE),
           list.files(file.path(root, "scripts"), pattern = "\\.[Rr]$", full.names = TRUE),
           list.files(file.path(root, "reference"), pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE))
result <- do.call(rbind, lapply(files, function(f) {
  error <- tryCatch({parse(f, encoding = "UTF-8"); ""}, error = function(e) conditionMessage(e))
  data.frame(file = substring(gsub("\\\\", "/", f), nchar(root) + 2L),
             syntax_ok = !nzchar(error), error = error)
}))
# Do not disclose the maintainer's absolute workstation paths in public reports.
result$error <- gsub(root, "<repo>", result$error, fixed = TRUE)
utils::write.csv(result, file.path(root, "docs/static_syntax_report.csv"), row.names = FALSE)
print(result[c("file", "syntax_ok")], row.names = FALSE)
if (any(!result$syntax_ok & !startsWith(result$file, "reference/")))
  stop("Formal module/entrypoint has a syntax error. No analysis was executed.")
cat("Static parse completed. Historical reference failures are listed, not silently fixed. No analysis executed.\n")
