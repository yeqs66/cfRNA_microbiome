# Run this AFTER the supplied script has created matrix_meta, tpm_1562,
# metabolics_matrix_val and select_pathway_micro. Do not source the original
# data-loading script from here. This file uses objects already in memory.
# Dependency (install yourself if needed): install.packages("igraph")
# Method: exploratory, unadjusted Spearman associations within each trimester
# plus an optional "all" group pooling every sample across trimesters.
# No imputation and no within-omics edges. The pooled analysis does NOT account
# for repeated participants or trimester effects; its p-values/FDR are nominal.

run_cross_omics_network <- function(
    matrix_meta, feature_sets,
    out_dir = file.path(getwd(), "matrix_meta_network",
                        format(Sys.time(), "%Y%m%d_%H%M%S")),
    min_n = 20L, max_missing = 0.30, min_micro_prevalence = 0.20,
    rho_cutoff = 0.30, fdr_cutoff = 0.05,
    p_adjust_method = "BH", max_plot_edges = 80L,
    label_top_n = 20L, seed = 2026L, include_all = TRUE) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop('Install igraph first: install.packages("igraph")')
  }
  stopifnot(min_n >= 3, max_missing >= 0, max_missing < 1,
            min_micro_prevalence >= 0, min_micro_prevalence <= 1,
            rho_cutoff >= 0, rho_cutoff <= 1,
            fdr_cutoff > 0, fdr_cutoff <= 1,
            max_plot_edges >= 1, label_top_n >= 0,
            p_adjust_method %in% stats::p.adjust.methods,
            is.logical(include_all), length(include_all) == 1L,
            !is.na(include_all))
  dat <- as.data.frame(matrix_meta, check.names = FALSE)
  if (!"sign" %in% names(dat) || anyDuplicated(names(dat))) {
    stop("matrix_meta must contain sign and have unique column names.")
  }
  dat$sign <- as.character(dat$sign)
  if (anyNA(dat$sign) || any(!nzchar(dat$sign)) || anyDuplicated(dat$sign)) {
    stop("Missing or duplicated sign. Check one-to-many joins upstream; do not arbitrarily drop duplicates.")
  }
  if (any(!grepl("^.+_(early|middle|late)$", dat$sign))) {
    stop("Expected sign format: ID_early, ID_middle or ID_late.")
  }
  # Remove only the FINAL suffix; IDs may themselves contain underscores.
  dat$ID_network <- sub("_(early|middle|late)$", "", dat$sign)
  dat$visit_network <- sub("^.*_(early|middle|late)$", "\\1", dat$sign)
  if (anyDuplicated(dat[c("ID_network", "visit_network")])) {
    stop("Each participant must have at most one row per trimester.")
  }

  omics <- c("Microbiome", "cfRNA", "Metabolome")
  if (!setequal(names(feature_sets), omics)) {
    stop("feature_sets must be named Microbiome, cfRNA and Metabolome.")
  }
  feature_sets <- lapply(feature_sets[omics], unique)
  if (any(lengths(feature_sets) == 0L)) stop("One omics feature set is empty.")
  all_features <- unlist(feature_sets, use.names = FALSE)
  if (anyNA(all_features) || any(!nzchar(all_features)) ||
      anyDuplicated(all_features)) {
    stop("Feature names must be unique across omics. Add omics prefixes BEFORE joining if names overlap.")
  }
  missing_features <- setdiff(all_features, names(dat))
  if (length(missing_features)) {
    stop("Features missing from matrix_meta (check .x/.y suffixes or renaming): ",
         paste(head(missing_features, 15), collapse = ", "))
  }
  if (any(all_features %in% c("sign", "ID_network", "visit_network"))) {
    stop("A feature name conflicts with a reserved metadata column.")
  }
  non_numeric <- all_features[!vapply(dat[all_features], is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop("Non-numeric features; inspect instead of coercing silently: ",
         paste(head(non_numeric, 15), collapse = ", "))
  }
  if (any(vapply(dat[feature_sets$Microbiome],
                 function(x) any(x < 0, na.rm = TRUE), logical(1)))) {
    stop("Microbiome prevalence assumes nonnegative RPK/abundance values.")
  }
  # NA, NaN and Inf are missing observations, never zero-filled.
  for (f in all_features) dat[[f]][!is.finite(dat[[f]])] <- NA_real_
  feature_map <- data.frame(
    feature = all_features, omics = rep(omics, lengths(feature_sets)),
    stringsAsFactors = FALSE)
  feature_map$node_id <- paste(feature_map$omics, feature_map$feature, sep = "::")
  node_id <- setNames(feature_map$node_id, feature_map$feature)
  visits <- c("early", "middle", "late")
  visits <- visits[visits %in% dat$visit_network]
  if (!length(visits)) stop("No usable samples.")
  if (include_all) visits <- c(visits, "all")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  write_csv <- function(x, filename) {
    utils::write.csv(x, file.path(out_dir, filename), row.names = FALSE,
                     na = "", fileEncoding = "UTF-8")
  }
  empty_edges <- function() data.frame(
    from = character(), to = character(), source_feature = character(),
    target_feature = character(), source_omics = character(),
    target_omics = character(), omics_pair = character(), visit = character(),
    n = integer(), n_participants = integer(),
    has_repeated_participants = logical(),
    rho = double(), p_value = double(),
    status = character(), stringsAsFactors = FALSE)
  edges_list <- list()
  qc_list <- list()
  sample_list <- list()
  omics_pairs <- combn(omics, 2, simplify = FALSE)

  for (v in visits) {
    # "all" uses the original sample rows exactly once; do not duplicate rows
    # by appending trimester-specific data before pooling.
    d <- if (v == "all") dat else dat[dat$visit_network == v, , drop = FALSE]
    n_people <- length(unique(d$ID_network))
    message("Processing ", v, ": ", nrow(d), " samples; ", n_people,
            " unique participants")
    if (v == "all" && n_people < nrow(d)) {
      message("Pooled Spearman treats sample rows as independent: repeated IDs are not modelled; p-values/FDR are nominal.")
    }
    has_omics <- sapply(feature_sets, function(fs)
      rowSums(!is.na(d[fs])) > 0)
    sample_list[[v]] <- data.frame(
      visit = v, n_samples = nrow(d), n_participants = n_people,
      has_repeated_participants = n_people < nrow(d),
      n_with_microbiome = sum(has_omics[, "Microbiome"]),
      n_with_cfRNA = sum(has_omics[, "cfRNA"]),
      n_with_metabolome = sum(has_omics[, "Metabolome"]),
      n_with_all_three = sum(rowSums(has_omics) == 3))

    # QC is recomputed within each analysis group, including pooled samples.
    # min_n and missingness thresholds count sample rows, not unique people.
    # Prevalence denominator is
    # the number of OBSERVED values for that feature, not all matrix rows.
    qc <- feature_map
    qc$visit <- v
    qc$n_observed <- vapply(qc$feature, function(f)
      sum(!is.na(d[[f]])), integer(1))
    qc$missing_fraction <- 1 - qc$n_observed / nrow(d)
    qc$n_unique <- vapply(qc$feature, function(f)
      length(unique(d[[f]][!is.na(d[[f]])])), integer(1))
    qc$prevalence <- vapply(qc$feature, function(f) {
      x <- d[[f]]
      if (!any(!is.na(x))) NA_real_ else mean(x > 0, na.rm = TRUE)
    }, numeric(1))
    qc$prevalence[qc$omics != "Microbiome"] <- NA_real_
    qc$keep <- qc$n_observed >= min_n &
      qc$missing_fraction <= max_missing & qc$n_unique >= 2 &
      (qc$omics != "Microbiome" |
         (!is.na(qc$prevalence) & qc$prevalence >= min_micro_prevalence))
    qc_list[[v]] <- qc

    for (pair in omics_pairs) {
      fs1 <- qc$feature[qc$keep & qc$omics == pair[1]]
      fs2 <- qc$feature[qc$keep & qc$omics == pair[2]]
      if (!length(fs1) || !length(fs2)) next
      pairs <- expand.grid(f1 = fs1, f2 = fs2, stringsAsFactors = FALSE)
      ans <- lapply(seq_len(nrow(pairs)), function(i) {
        f1 <- pairs$f1[i]; f2 <- pairs$f2[i]
        ok <- !is.na(d[[f1]]) & !is.na(d[[f2]])
        x <- d[[f1]][ok]; y <- d[[f2]][ok]
        n <- length(x)
        n_participants <- length(unique(d$ID_network[ok]))
        rho <- p <- NA_real_
        status <- "ok"
        if (n < min_n) {
          status <- "insufficient_pairs"
        } else if (length(unique(x)) < 2 || length(unique(y)) < 2) {
          status <- "constant_in_paired_samples"
        } else {
          # exact=FALSE handles tied ranks; p-values are asymptotic.
          ct <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
          rho <- unname(ct$estimate); p <- ct$p.value
          if (!is.finite(rho) || !is.finite(p)) status <- "nonfinite_test"
        }
        data.frame(
          from = unname(node_id[f1]), to = unname(node_id[f2]),
          source_feature = f1, target_feature = f2,
          source_omics = pair[1], target_omics = pair[2],
          omics_pair = paste(pair, collapse = "--"), visit = v,
          n = n, n_participants = n_participants,
          has_repeated_participants = n_participants < n,
          rho = rho, p_value = p, status = status,
          stringsAsFactors = FALSE)
      })
      edges_list[[paste(v, paste(pair, collapse = "_"), sep = "_")]] <-
        do.call(rbind, ans)
    }
  }

  edges <- if (length(edges_list)) do.call(rbind, edges_list) else empty_edges()
  rownames(edges) <- NULL
  edges$FDR <- NA_real_
  valid <- edges$status == "ok" & is.finite(edges$p_value)
  # Keep the previous testing family for early/middle/late together. The
  # additional pooled analysis has its own family across all omics pairs.
  # Adding "all" therefore does not change the trimester-specific FDR values.
  # Pooled FDR is nominal if sample rows include repeated participants.
  # Adjust BEFORE applying the effect-size threshold or selecting plot edges.
  edges$FDR_family <- ifelse(edges$visit == "all", "pooled_all", "trimesters_joint")
  for (family in unique(edges$FDR_family)) {
    idx <- valid & edges$FDR_family == family
    edges$FDR[idx] <- stats::p.adjust(edges$p_value[idx], method = p_adjust_method)
  }
  edges$analysis_note <- ifelse(edges$visit == "all",
    "Pooled unadjusted Spearman; no repeated-measures or trimester adjustment; nominal p/FDR",
    "Within-trimester unadjusted Spearman")
  edges$abs_rho <- abs(edges$rho)
  edges$weight <- edges$abs_rho  # positive layout weights, NOT signed rho
  edges$direction <- ifelse(is.na(edges$rho), NA_character_,
                            ifelse(edges$rho >= 0, "positive", "negative"))
  significant <- valid & !is.na(edges$FDR) & edges$FDR < fdr_cutoff &
    edges$abs_rho >= rho_cutoff
  sig <- edges[significant, , drop = FALSE]
  sig <- sig[order(sig$FDR, -sig$abs_rho, sig$from, sig$to), , drop = FALSE]
  qc_all <- do.call(rbind, qc_list)
  write_csv(edges, "associations_all.csv")
  write_csv(sig, "edges_significant_all_visits.csv")
  write_csv(qc_all, "feature_QC.csv")
  write_csv(feature_map, "feature_map.csv")
  write_csv(do.call(rbind, sample_list), "sample_summary.csv")

  node_colors <- c(Microbiome = "#009E73", cfRNA = "#0072B2",
                   Metabolome = "#E69F00")
  node_shapes <- c(Microbiome = "square", cfRNA = "circle",
                   Metabolome = "csquare")
  edge_colors <- c(positive = "#B35806", negative = "#542788")
  vertices_for <- function(e) {
    ids <- unique(c(e$from, e$to))
    nd <- feature_map[match(ids, feature_map$node_id), , drop = FALSE]
    data.frame(name = nd$node_id, label = nd$feature, omics = nd$omics,
               color = unname(node_colors[nd$omics]),
               shape = unname(node_shapes[nd$omics]), stringsAsFactors = FALSE)
  }
  # Select only the plotting subset; full significant networks remain exported.
  plot_lists <- lapply(visits, function(v) {
    ev <- sig[sig$visit == v, , drop = FALSE]
    head(ev, max_plot_edges)
  })
  names(plot_lists) <- visits
  plot_edges <- do.call(rbind, plot_lists)
  write_csv(plot_edges, "edges_displayed_all_visits.csv")

  # Use one common layout for trimester-specific and pooled figures.
  if (nrow(plot_edges)) {
    union_edges <- unique(plot_edges[c("from", "to")])
    union_graph <- igraph::graph_from_data_frame(
      union_edges, directed = FALSE, vertices = vertices_for(union_edges))
    set.seed(seed)
    common_xy <- igraph::layout_with_fr(union_graph)
    common_xy <- igraph::norm_coords(common_xy, xmin = -1, xmax = 1,
                                     ymin = -1, ymax = 1)
    rownames(common_xy) <- igraph::V(union_graph)$name
  }

  node_lists <- list()
  summaries <- list()
  for (v in visits) {
    ev <- sig[sig$visit == v, , drop = FALSE]
    pe <- plot_lists[[v]]
    nd <- vertices_for(ev)
    nd$visit <- rep(v, nrow(nd))
    if (nrow(ev)) {
      g <- igraph::graph_from_data_frame(ev, directed = FALSE, vertices = nd)
      nd$degree <- as.numeric(igraph::degree(g)[nd$name])
      nd$strength <- as.numeric(igraph::strength(g, weights = igraph::E(g)$weight)[nd$name])
      igraph::V(g)$degree <- nd$degree
      igraph::V(g)$strength <- nd$strength
      igraph::write_graph(g, file.path(out_dir, paste0("network_", v, ".graphml")),
                          format = "graphml")
    } else {
      nd$degree <- integer(0); nd$strength <- numeric(0)
    }
    node_lists[[v]] <- nd
    write_csv(ev, paste0("edges_", v, ".csv"))
    write_csv(nd, paste0("nodes_", v, ".csv"))
    write_csv(pe, paste0("edges_displayed_", v, ".csv"))
    summaries[[v]] <- data.frame(
      visit = v, n_tested = sum(valid & edges$visit == v),
      n_significant_edges = nrow(ev), n_network_nodes = nrow(nd),
      n_displayed_edges = nrow(pe),
      FDR_family = if (v == "all") "pooled_all" else "trimesters_joint")

    draw <- function() {
      graphics::par(mar = c(3.5, 2, 5, 2), family = "sans")
      plot_title <- if (v == "all") "Cross-omics network: all trimesters pooled" else
        paste("Cross-omics network:", v)
      pooled_note <- "Unadjusted pooled samples; repeated participants not modelled; nominal p/FDR."
      if (!nrow(pe)) {
        graphics::plot.new()
        graphics::title(main = plot_title)
        graphics::text(0.5, 0.5, "No associations pass the prespecified thresholds.")
        if (v == "all") graphics::mtext(pooled_note, side = 1, line = 1.5, cex = 0.75)
        return(invisible(NULL))
      }
      gp <- igraph::graph_from_data_frame(pe, directed = FALSE,
                                         vertices = vertices_for(pe))
      deg <- igraph::degree(gp)
      label_ids <- names(sort(deg, decreasing = TRUE))[seq_len(min(label_top_n, length(deg)))]
      # Label selected hubs and metabolites; complete labels remain in node tables.
      labels <- ifelse(igraph::V(gp)$name %in% label_ids |
                         igraph::V(gp)$omics == "Metabolome",
                       igraph::V(gp)$label, "")
      labels <- vapply(labels, function(z)
        paste(strwrap(z, width = 26), collapse = "\n"), character(1))
      xy <- common_xy[igraph::V(gp)$name, , drop = FALSE]
      plot(gp, layout = xy, rescale = FALSE,
           xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45),
           vertex.size = 5 + 2 * sqrt(deg), vertex.label = labels,
           vertex.label.cex = 0.65, vertex.label.color = "#222222",
           vertex.label.dist = 0.7, vertex.frame.color = "white",
           edge.color = grDevices::adjustcolor(
             unname(edge_colors[igraph::E(gp)$direction]), alpha.f = 0.65),
           edge.lty = ifelse(igraph::E(gp)$direction == "positive", 1, 2),
           edge.width = 0.6 + 3 * igraph::E(gp)$abs_rho,
           edge.curved = 0.06, edge.arrow.mode = 0)
      graphics::title(main = plot_title,
        sub = sprintf("%s FDR < %.3g; |rho| >= %.2f; %d/%d significant edges displayed",
                      p_adjust_method, fdr_cutoff, rho_cutoff, nrow(pe), nrow(ev)),
        cex.sub = 0.8)
      if (v == "all") graphics::mtext(pooled_note, side = 1, line = 2.5, cex = 0.75)
      graphics::legend("topleft", legend = omics, pch = c(15, 16, 0),
                       col = node_colors[omics], bty = "n", cex = 0.85)
      graphics::legend("topright", legend = c("Positive", "Negative"),
                       col = edge_colors, lty = c(1, 2), lwd = 2,
                       bty = "n", cex = 0.85)
      invisible(NULL)
    }
    save_plot <- function(device, filename, ...) {
      device(filename, ...)
      on.exit(grDevices::dev.off(), add = TRUE)
      draw()
    }
    save_plot(grDevices::pdf, file.path(out_dir, paste0("network_", v, ".pdf")),
              width = 14, height = 11, useDingbats = FALSE)
    save_plot(grDevices::png, file.path(out_dir, paste0("network_", v, ".png")),
              width = 4200, height = 3300, res = 300)
  }
  nodes <- do.call(rbind, node_lists)
  summary_table <- do.call(rbind, summaries)
  write_csv(nodes, "nodes_significant_all_visits.csv")
  write_csv(summary_table, "network_summary.csv")
  settings <- list(method = "Spearman; asymptotic two-sided p-values",
    strata = visits, include_all = include_all,
    adjustment = "unadjusted for clinical covariates",
    pooled_analysis = "all original sample rows pooled; no repeated-measures or trimester adjustment; nominal p/FDR",
    multiple_testing_family = c(
      trimesters_joint = "all valid cross-omics tests across early/middle/late jointly",
      pooled_all = "all valid cross-omics tests in the pooled analysis, separately"),
    p_adjust_method = p_adjust_method, min_n = min_n,
    max_missing = max_missing, min_micro_prevalence = min_micro_prevalence,
    rho_cutoff = rho_cutoff, fdr_cutoff = fdr_cutoff,
    display_selection = "lowest FDR, then largest absolute rho, within each visit",
    max_plot_edges = max_plot_edges, seed = seed)
  writeLines(capture.output(dput(settings)), file.path(out_dir, "analysis_settings.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  result <- list(associations = edges, edges = sig, nodes = nodes,
                 feature_QC = qc_all, summary = summary_table,
                 settings = settings, out_dir = out_dir)
  saveRDS(result, file.path(out_dir, "network_results.rds"))
  message("Output directory: ", out_dir)
  print(summary_table)
  invisible(result)
}

# ------------------ Run after building matrix_meta ------------------
# Define features from the ORIGINAL omics tables; avoid classifying every
# unrecognised column as a metabolite (clinical variables could be included).
feature_sets <- list(
  Microbiome = unique(as.character(select_pathway_micro)),
  cfRNA = setdiff(colnames(tpm_1562), "sign"),
  Metabolome = setdiff(colnames(metabolics_matrix_val), "sign")
)
# select_pathway_micro may contain pathway KOs absent from the source table;
# analyse only KOs actually present in matrix_meta, matching the original filter.
feature_sets$Microbiome <- intersect(feature_sets$Microbiome, colnames(matrix_meta))

matrix_network<-matrix_meta_log %>% dplyr::filter(sign %in% coldata_all$sign)

network_result <- run_cross_omics_network(
  matrix_meta = matrix_network,
  feature_sets = feature_sets,
  min_n = 20,
  max_missing = 0.50,
  min_micro_prevalence = 0.10,
  rho_cutoff = 0.10,
  fdr_cutoff = 0.05,
  p_adjust_method = "BH",  # use "BY" for arbitrary-dependence FDR control
  max_plot_edges = 80,
  label_top_n = 20,
  include_all = TRUE  # early, middle, late AND all samples pooled
)

# View results in R:
# network_result$summary
# head(network_result$edges)
# network_result$out_dir
# subset(network_result$edges, visit == "all")
# Set include_all = FALSE to run only trimester-specific analyses.
# To display all significant edges: set max_plot_edges = Inf above.
# For GDM/control-specific analyses, first subset matrix_meta by signs taken
# from the appropriate group in coldata_all, then rerun into a NEW output folder.
