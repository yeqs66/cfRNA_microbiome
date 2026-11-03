# Extract the ORIGINAL fitting functions without evaluating data-loading/plot blocks.
run_metabolite_lmm <- function(data, metabolite_map, out_dir, repo_root) {
  require_packages(c("dplyr", "lme4", "lmerTest", "emmeans"))
  validate_samples(data)
  covariates <- c("batch", "BMI_mo_Anthropometric", "ethnicity_mo",
                 "smoking_mo", "drinking_mo_18topreg")
  required <- c("sign", "ID", "trimester", "disease", covariates)
  if (length(setdiff(required, names(data)))) stop("Missing metadata: ",
    paste(setdiff(required, names(data)), collapse = ", "))
  if (!all(c("Compound", "category") %in% names(metabolite_map)) ||
      anyNA(metabolite_map$Compound) || anyDuplicated(metabolite_map$Compound) ||
      !nrow(metabolite_map)) stop("Provide a nonempty unique Compound/category map.")
  features <- as.character(metabolite_map$Compound)
  if (length(setdiff(features, names(data)))) stop("Metabolite columns absent.")
  if (!all(vapply(data[features], is.numeric, logical(1)))) stop("Metabolites must be numeric.")
  if (any(vapply(data[features], function(x) any(is.infinite(x)), logical(1))))
    stop("Infinite concentrations: inspect input before fitting.")
  if (anyNA(data$ID) || any(!nzchar(as.character(data$ID)))) stop("ID missing.")
  if (any(!as.character(data$disease) %in% c("GDM", "Control")) || anyNA(data$disease))
    stop("disease must be GDM/Control.")
  if (anyNA(data$trimester) || any(!as.character(data$trimester) %in% c("early", "middle", "late")))
    stop("Invalid trimester.")
  expected_sign <- paste(data$ID, data$trimester, sep = "_")
  if (!identical(as.character(data$sign), expected_sign)) stop("sign/ID/trimester mismatch.")
  data$disease <- factor(data$disease, levels = c("GDM", "Control"))
  data$trimester <- factor(data$trimester, levels = c("early", "middle", "late"))
  for (f in setdiff(covariates, "BMI_mo_Anthropometric")) data[[f]] <- factor(data[[f]])
  env <- source_functions(file.path(repo_root, "reference/metabolomics/01.R"),
    c("log2fc_fun", "form_primary", "form_adjusted", "fit_one"))
  # Namespace the original function calls; formulas and contrasts are unchanged.
  env$`%>%` <- dplyr::`%>%`
  env$transmute <- dplyr::transmute
  env$bind_rows <- dplyr::bind_rows
  env$lmer <- lmerTest::lmer
  env$emmeans <- emmeans::emmeans
  env$metabolics_matrix_val <- data
  env$adjusted_cov <- covariates
  result <- dplyr::bind_rows(lapply(c("primary", "adjusted"), function(m)
    dplyr::bind_rows(lapply(features, env$fit_one, model = m))))
  result <- dplyr::group_by(result, model, trimester)
  result <- dplyr::mutate(result, padj = p.adjust(pvalue, method = "fdr"),
                           logp = -log10(padj))
  result <- dplyr::ungroup(result)
  result <- dplyr::left_join(result, metabolite_map[c("Compound", "category")], by = "Compound")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(out_dir, "metabolite_lmm.csv"), row.names = FALSE, na = "")
  utils::write.csv(as.data.frame(table(data$trimester, data$disease)),
    file.path(out_dir, "sample_summary.csv"), row.names = FALSE)
  # Original 'n' counts the full fitted model, not the number in a contrast's trimester.
  saveRDS(list(primary = env$form_primary(features[1]),
    adjusted = env$form_adjusted(features[1]), covariates = covariates,
    all_contrast = "equal-weight marginal contrast across three trimesters",
    FDR = "BH within model x trimester across supplied metabolite map",
    scale = "untransformed concentrations; descriptive log2FC is not model estimate",
    features = features), file.path(out_dir, "analysis_settings.rds"))
  saveRDS(result, file.path(out_dir, "results.rds"))
  result
}
