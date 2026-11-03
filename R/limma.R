# Refactor of reference/cfrna/01.R and 02.R; retain original contrasts/filtering.
# See audit: these contrasts and the >1 filtering differ from manuscript wording.
run_limma_reference <- function(counts, metadata, rna_type, out_dir) {
  require_packages(c("limma", "edgeR"))
  stopifnot(rna_type %in% c("mirna", "lncrna"))
  validate_samples(metadata)
  required <- c("sample", "ID", "trimester", "disease", "batch")
  if (length(setdiff(required, names(metadata)))) stop("Missing limma metadata.")
  if (anyDuplicated(metadata$sample) || anyNA(metadata$sample) || anyNA(metadata$ID))
    stop("Sample/ID missing or sample duplicated.")
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      is.null(rownames(counts)) || anyDuplicated(rownames(counts)) ||
      is.null(colnames(counts)) || anyDuplicated(colnames(counts))) stop("Invalid counts matrix.")
  if (length(setdiff(metadata$sample, colnames(counts)))) stop("Counts missing requested samples.")
  counts <- counts[, metadata$sample, drop = FALSE]
  # Exact historical condition; do not silently change to >=1 in >=10%.
  keep <- rowSums(counts > 1) > floor(ncol(counts) / 10)
  counts <- counts[keep, , drop = FALSE]
  if (nrow(counts) < 2L) stop("Insufficient genes after original filtering.")
  metadata$trimester <- factor(as.character(metadata$trimester), levels = c("early", "middle", "late"))
  metadata$disease <- factor(as.character(metadata$disease), levels = c("Control", "GDM"))
  metadata$batch <- factor(metadata$batch)
  if (anyNA(metadata[required])) stop("Missing/invalid limma metadata values.")
  design <- model.matrix(~trimester * disease + batch, data = metadata)
  if (ncol(design) != 7L || qr(design)$rank != 7L)
    stop("Historical contrasts require seven full-rank coefficients (three trimesters, two diseases, two batches). See audit.")
  contrasts <- cbind(de_1 = c(0,0,0,1,0,0,0), de_2 = c(0,0,0,1,0,1,0),
    de_3 = c(0,0,0,1,0,1,1), de_all = c(0,0,0,1,0,1/3,1/3))
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
  w0 <- limma::voomWithQualityWeights(dge, design, normalize.method = "none",
                                     method = "genebygene", trace = TRUE)
  corfit <- limma::duplicateCorrelation(w0, design, block = metadata$ID)
  # Original second voom does not pass block. Preserve and flag in audit.
  w <- limma::voomWithQualityWeights(dge, design, normalize.method = "none",
    method = "genebygene", correlation = corfit$consensus, trace = TRUE)
  fit <- limma::lmFit(w, design, block = metadata$ID, correlation = corfit$consensus)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result <- lapply(colnames(contrasts), function(coef) {
    tab <- limma::topTable(fit, coef = coef, adjust.method = "BH", n = Inf)
    tab$feature <- rownames(tab)
    tab$period <- sub("de_", "", coef)
    utils::write.csv(tab, file.path(out_dir, paste0(coef, ".csv")), row.names = FALSE)
    tab
  })
  names(result) <- colnames(contrasts)
  utils::write.csv(data.frame(feature = rownames(as.matrix(counts))),
    file.path(out_dir, "retained_features.csv"), row.names = FALSE)
  saveRDS(list(design = design, contrasts = contrasts, consensus = corfit$consensus,
               filter = ">1 counts in >floor(n_samples/10) samples", rna_type = rna_type),
          file.path(out_dir, "analysis_settings.rds"))
  saveRDS(result, file.path(out_dir, "results.rds"))
  result
}
