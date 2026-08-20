# Data processing + anonymization (mock pipeline)
library(dplyr)
library(tidyr)
library(MatrixGenerics)
library(sva)

# ------------------------------ 1. Simulate raw data ------------------------------
set.seed(42)
n_genes <- 500
n_samples <- 100
raw_counts <- matrix(rnbinom(n_genes * n_samples, mu = 10, size = 2),
                     nrow = n_genes, ncol = n_samples,
                     dimnames = list(paste0("gene", 1:n_genes),
                                     paste0("sample", 1:n_samples)))

# clinical metadata
coldata <- data.frame(
  sample = colnames(raw_counts),
  disease = sample(c(0,1), n_samples, replace = TRUE),
  batch = sample(c("batch1","batch2"), n_samples, replace = TRUE),
  trimester = sample(c("1","2","3"), n_samples, replace = TRUE),
  ID = paste0("ID", sample(100:199, n_samples, replace = TRUE)),
  stringsAsFactors = FALSE
)

# microbiome KO matrix (mock)
ko_counts <- matrix(rnbinom(200 * n_samples, mu = 5, size = 1.5),
                    nrow = 200, ncol = n_samples,
                    dimnames = list(paste0("KO", 1:200), colnames(raw_counts)))

# ------------------------------ 2. Filter low-expressed genes ------------------------------
filter_low <- function(mat, min_samps = 10, min_count = 2) {
  keep <- rowSums(mat > min_count) >= min_samps
  mat[keep, , drop = FALSE]
}
rna_flt <- filter_low(raw_counts, min_samps = floor(ncol(raw_counts)/2))
ko_flt <- filter_low(ko_counts, min_samps = floor(ncol(ko_counts)/10))

# ------------------------------ 3. Variance filter (top 50%) ------------------------------
filter_top_var <- function(mat, top_frac = 0.5) {
  vars <- rowVars(mat)
  thresh <- quantile(vars, 1 - top_frac, na.rm = TRUE)
  mat[vars >= thresh, , drop = FALSE]
}
rna_var <- filter_top_var(rna_flt, top_frac = 0.5)
ko_var <- filter_top_var(ko_flt, top_frac = 0.5)

# ------------------------------ 4. Log2 transform (add pseudocount) ------------------------------
rna_log <- log2(rna_var + 1)
ko_log <- log2(ko_var + 1)

# ------------------------------ 5. ComBat batch correction (if batch present) ------------------------------
if ("batch" %in% colnames(coldata)) {
  batch <- as.factor(coldata$batch)
  mod <- model.matrix(~ disease, data = coldata)  # preserve disease effect
  rna_combat <- ComBat(rna_log, batch = batch, mod = mod, par.prior = TRUE)
  ko_combat <- ComBat(ko_log, batch = batch, mod = mod, par.prior = TRUE)
} else {
  rna_combat <- rna_log
  ko_combat <- ko_log
}

# align samples (ensure same order)
common_samples <- intersect(colnames(rna_combat), colnames(ko_combat))
coldata <- coldata[match(common_samples, coldata$sample), , drop = FALSE]
rna_combat <- rna_combat[, common_samples, drop = FALSE]
ko_combat <- ko_combat[, common_samples, drop = FALSE]

# ------------------------------ 6. Build MOFA input list ------------------------------
mofa_input <- list(
  mRNA = rna_combat,
  KO   = ko_combat
)
# Add other views if needed (miRNA, lncRNA) – mock
# mofa_input$miRNA <- ...
# mofa_input$lncRNA <- ...

# ------------------------------ 7. Anonymization functions ------------------------------
anonymize <- function(ids, prefix = "SUBJ", seed = 123) {
  set.seed(seed)
  n <- length(ids)
  perm <- sample(n)
  new_ids <- paste0(prefix, "_", sprintf(paste0("%0", nchar(n), "d"), perm))
  return(new_ids)
}

# generate mapping
id_map <- data.frame(
  original = coldata$ID,
  anon_id = anonymize(coldata$ID)
)

# apply to sample names (if needed)
anon_sample <- setNames(id_map$anon_id, coldata$sample)
colnames(rna_combat) <- anon_sample[colnames(rna_combat)]
colnames(ko_combat) <- anon_sample[colnames(ko_combat)]

# update coldata with anonymized IDs
coldata_anon <- coldata
coldata_anon$sample <- anon_sample[coldata$sample]
coldata_anon$ID <- id_map$anon_id[match(coldata$ID, id_map$original)]

# ------------------------------ 8. Output (mock save) ------------------------------
# In real use, write tables to file:
# write.csv(rna_combat, "rna_processed.csv")
# write.csv(ko_combat, "ko_processed.csv")
# write.csv(coldata_anon, "coldata_anon.csv")
# write.csv(id_map, "id_mapping.csv")

# Show summaries
cat("Processed RNA matrix:", dim(rna_combat), "\n")
cat("Processed KO matrix:", dim(ko_combat), "\n")
cat("MOFA input list views:", names(mofa_input), "\n")
cat("Anonymization mapping (first 5):\n")
print(head(id_map))