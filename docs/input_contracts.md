# Inputs and execution contracts

No study inputs are included or opened during organization. Configure private files yourself; relative paths are resolved against the YAML configuration directory. All outputs are written to a new timestamp/PID/stage directory. No analysis runs without explicit `--stage` and `--mode`.

## Common sample keys

- `sign`: unique `ID_early`, `ID_middle`, `ID_late`; IDs may contain underscores. Duplicate keys are errors, never silently dropped.
- Covariates/outcomes must be matched by key, not row order. A single participant can have one record in each trimester.
- Network features must be numeric with unique names across omics layers. Only the explicit feature lists are tested; metadata must not be classified as metabolites.

## network

Input CSV or RDS data.frame: rows are samples; columns are `sign` plus numeric features. `feature_sets` names must be Microbiome/cfRNA/Metabolome. For the manuscript use the author-confirmed nine metabolites, gene annotations and KO membership, not the illustrative two-feature lists in the configuration. Microbiome values must be nonnegative. Missing/Inf values are treated by the preserved source function; no imputation. Log transformation of nonnegative KO abundances preserves nonzero prevalence but must match the final study pipeline.

Output: all tested edges, valid sample/participant counts, feature QC, significant edges, separately selected display edges, nodes, GraphML, PDF/PNG, analysis settings and RDS. Early/middle/late use joint BH; pooled `all` uses separate BH and does not adjust repeated observations. `max_plot_edges` affects displays only.

## metabolomics_lmm

Input CSV/RDS data.frame: `sign`, `ID`, `trimester`, `disease` (Control/GDM), `batch`, `BMI_mo_Anthropometric` (numeric), `ethnicity_mo`, `smoking_mo`, `drinking_mo_18topreg`, and numeric Compound columns. Values are untransformed concentrations, as in source. Metabolite map has unique `Compound`, `category` columns. Its membership determines the FDR family; do not restrict it to nine metabolites without resolving audit M01.

Output includes both primary and adjusted models, period-specific and equal-weight all-period contrasts, original descriptive log2FC and BH within model×trimester. Original fitting failures produce NA rows/messages; singular-fit warnings were suppressed by source, see audit. `n` is full model fitting sample count. This module deliberately does not promise the manuscript's exact results.

## limma_mirna / limma_lncrna

Counts RDS matrix: rows unique genes, columns unique `sample` names, finite nonnegative numeric counts. Metadata CSV/RDS: `sign`, `sample`, `ID`, `trimester` early/middle/late, `disease` Control/GDM, `batch`. Counts are reordered to metadata sample names. Historical filtering and fixed contrasts are preserved; metadata must yield a full-rank seven-column design. In real mode `acknowledge_legacy_limma: true` is required after reviewing audit L02–L04. This is not a substitute for the missing confirmed mRNA primary pipeline.

## Other manuscript stages

44 curated source candidates are organized under reference/; paths map back to the original scripts using docs/script_inventory.json. These are not standalone validated pipelines. Use the opt-in reviewed reference adapter only after completing its inputs, path mappings and upstream objects. A `globals_rds` file is a named list of precomputed objects; the adapter never reads an unconfigured private global environment.

`config/reference_path_tokens.json` lists per-script placeholders and suggested file basenames (not original private absolute paths). Configure token keys under a stage's `paths`. Map remaining relative input literals explicitly when they are not relative to the new output directory. The adapter preflights all redacted path tokens before evaluation, disables setwd/workspace-clearing calls and incidental save.image dumps, and evaluates the remaining historical expressions in order in an isolated environment. It does not resolve missing symbols, CLI positional arguments, reference resources, data units or analysis choices for you. Scripts needing these remain pending manual review.

No raw sequencing, deconvolution, GTDB-Tk or online TissueEnrich execution implementation is fabricated. Their upstream outputs may be supplied with documented versions after author review.
