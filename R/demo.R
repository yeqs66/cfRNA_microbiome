make_demo <- function() {
  set.seed(2026)
  n <- 36L
  id <- rep(sprintf("demo_%03d", seq_len(n)), 3)
  trimester <- rep(c("early", "middle", "late"), each = n)
  disease <- rep(rep(c("Control", "GDM"), each = n / 2), 3)
  subject <- rep(rnorm(n), 3)
  z <- rnorm(length(id))
  network <- data.frame(sign = paste(id, trimester, sep = "_"),
    KO_001 = exp(z), KO_002 = exp(-z + rnorm(length(id), sd = .2)),
    RNA_001 = z + rnorm(length(id), sd = .2),
    RNA_002 = -z + rnorm(length(id), sd = .2), check.names = FALSE)
  metabolites <- c("L-Valine(Val)", "Isoleucine", "L-Leucine(Leu)",
    "alpha-ketoisovaleric acid", "ketoleucine", "3-methyl-2-oxovaleric acid",
    "3-hydroxyisovaleric acid", "acetoacetic acid", "methylmalonic acid")
  metadata <- data.frame(sign = network$sign, ID = id,
    sample = paste0("demoRNA_", seq_along(id)), trimester = trimester,
    disease = disease, batch = rep(c("batchA", "batchB"), length.out = length(id)),
    BMI_mo_Anthropometric = rep(runif(n, 18, 30), 3),
    ethnicity_mo = rep(sample(c("A", "B"), n, TRUE), 3),
    smoking_mo = rep(sample(c("no", "yes"), n, TRUE), 3),
    drinking_mo_18topreg = rep(sample(c("no", "yes"), n, TRUE), 3))
  lmm <- metadata
  for (i in seq_along(metabolites)) {
    network[[metabolites[i]]] <- z + rnorm(length(id), sd = .4 + i / 10)
    lmm[[metabolites[i]]] <- 50 + 3 * (disease == "GDM") + subject +
      as.integer(factor(trimester)) + rnorm(length(id))
  }
  network$RNA_002[1:2] <- NA_real_
  counts <- matrix(rnbinom(60 * length(id), size = 10, mu = 100), nrow = 60,
    dimnames = list(paste0("demo_feature_", 1:60), metadata$sample))
  list(network = network, lmm = lmm, counts = counts, metadata = metadata,
       metabolite_map = data.frame(Compound = metabolites, category = "BCAA"),
       feature_sets = list(Microbiome = c("KO_001", "KO_002"),
                           cfRNA = c("RNA_001", "RNA_002"), Metabolome = metabolites))
}
