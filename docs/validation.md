# Validation record

## Execution boundary — 2026-11-03

At the author's request: no reproduction, study-data loading, statistical model execution or simulation analysis. Only source/document reading, static parsing and Git checks are performed.

An earlier network-entrypoint attempt before this restriction failed when sourcing infrastructure because of Windows locale/path handling, before statistical analysis executed. Locale handling has since been edited; the entrypoint is not rerun. Synthetic fixtures and future run examples are supplied without claiming dynamic validation.

## Checks performed

- Base R parse() on six modules, two entry/check scripts and 44 source references: 52 syntax successes; see static_syntax_report.csv.
- 139 source scripts inventoried by path/hash; 54 damaged/undecodable excluded; 44 selected source references retained after private-path redaction.
- Original main backed up as archive/pre-v9-20261102 at expected commit 4cf58be1f88b5eb21b878c8a55da1589003ef4f8; remote equality verified before publication.
- Original local script checksums verified against inventory before publication; no originals intentionally edited.
- Repository additions inspected to exclude data matrices, manuscript files and identity maps.

## Not validated

Runtime correctness, numerical equivalence, package compatibility, real input schemas, model convergence, memory/GPU needs and manuscript numerical claims remain untested. Source recipes and missing methods require author review. Static parsing alone cannot establish reproducibility.
