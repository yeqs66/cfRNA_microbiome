# cfRNA–microbiome: manuscript v9 code organization

**Longitudinal Multi-omics Profiling Reveals Systematic Host–Microbiome Metabolic Dysregulation in Gestational Diabetes Mellitus**

This is a source-backed, **partially configured** manuscript research project. At the author's request, validation is **static only**: no study data are opened and no models, simulations or numerical reproduction are run. Do not batch-run historical scripts on a personal workstation.

请先阅读：[论文与代码不一致核对清单（中文）](docs/manuscript_v9_code_audit.md)。缺失或不一致的实现留待作者确认，没有按稿件补写算法或调整统计结果。

## Contents and provenance

- [Inventory](docs/script_inventory.md) and [detailed records](docs/script_inventory.json): 139 source R scripts, original SHA-256, input/output references, packages and selection status.
- 54 damaged/undecodable sources excluded; 44 selected candidates in `reference/`, with private absolute paths replaced by configuration placeholders. These are review materials, not confirmed final versions for every figure.
- Independent interfaces for the original network function, original metabolite LMM functions, and historical miRNA/lncRNA limma code. Legacy limma contrasts and filtering require review before real-data use.
- An opt-in adapter for reviewed historical sources with configured paths/upstream objects. Missing sequencing, deconvolution and nested-CV implementations are not fabricated.
- Old repository code retained unchanged under `legacy/pre_v9/`. The complete old main is preserved on `archive/pre-v9-20261102`, commit `4cf58be1f88b5eb21b878c8a55da1589003ef4f8`.

## Static review

```sh
Rscript scripts/check_syntax.R
```

This parses code without evaluating it or loading study matrices. See [static report](docs/static_syntax_report.csv). The pipeline has no implicit run-all command and does not run without explicit stage/mode arguments.

## Future optional execution — not tested in this release

Read [input contracts](docs/input_contracts.md), [coverage](docs/analysis_coverage.md), and the audit before execution. Copy `config/config.example.yml` to ignored `config/local.yml`, configure authorized data paths and author-confirmed feature membership.

```sh
# Future small synthetic example; NOT executed as part of this delivery
Rscript scripts/run_pipeline.R --config config/demo.yml --stage network --mode demo

# Future real-data run after audit and data review
Rscript scripts/run_pipeline.R --config config/local.yml --stage network --mode real
Rscript scripts/run_pipeline.R --config config/local.yml --stage metabolomics_lmm --mode real
```

Other standalone stages: `limma_mirna`, `limma_lncrna`. Real limma runs require `acknowledge_legacy_limma: true`. Unknown/missing stages stop clearly. Reference recipes require `reviewed: true`, paths and upstream objects; that flag is not a claim of reproducibility.

Each future run uses a separate output directory and records mode, configuration, completion/failure and software environment. The network exports all tests/significant edges separately from its display subset. Pooled correlations do not adjust repeated observations.

## Dependencies, data and limitations

See [dependencies](docs/dependencies.md) and [validation record](docs/validation.md). Packages/resources are not installed or downloaded automatically. No guessed historical `renv.lock` is supplied. Original data, private clinical records, identities, manuscript DOCX and pretrained models are not bundled.

The new wrappers have passed **syntax parsing only**. Runtime behavior, wrapper/source numerical equivalence, resource needs and manuscript AUC/R²/P/FDR remain unverified. `reference/` line order is retained; original paths/hashes map through the inventory. Source candidates may contain obsolete blocks requiring manual review.

The previous README's blanket description of all code as simulated does not apply to the new project. This organization does not change licensing: the previous README mentioned MIT but the inspected tree did not contain a LICENSE file; confirm licensing before further redistribution.
