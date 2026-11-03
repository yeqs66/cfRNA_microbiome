# Dependencies and environment

No installation is performed automatically. Manuscript versions are targets for review, not a reconstructed lockfile.

| Stage | Dependencies |
|---|---|
| Static syntax check | Base R only; parsing checked on R4.5.2 |
| YAML runner | yaml |
| Network | igraph, base R |
| Metabolite LMM | dplyr, lme4, lmerTest, emmeans |
| Historical limma | edgeR, limma; statmod may be required by underlying routines |
| WGCNA/preparation | WGCNA, DESeq2, limma/sva, tidyverse; source-specific imports |
| Enrichment | clusterProfiler, graphite, org.Hs.eg.db, GTF/GMT; historical KEGG.db usage needs review |
| MOFA | MOFA2, reticulate, compatible Python mofapy2/basilisk and HDF5 |
| miRNA targets | multiMiR and confirmed database access/version |
| Prediction references | keras/tensorflow, compatible Python TensorFlow, pROC; other imports in inventory |
| Figures | ggplot2, ComplexHeatmap, circlize, patchwork, openxlsx/readxl; per-source imports |

The manuscript reports edgeR3.38.4, limma3.52.4, clusterProfiler4.9.0.002, WGCNA1.72-1 and glmmTMB1.1.9 among its dependencies. They have not been installed or validated as a compatible historical set. Review original session logs/database versions before pinning an environment. Current package availability is not evidence that a version generated the original results.

Some imports are unused or from alternative branches. Install only the reviewed stage you intend to run. No model execution, dependency installation or reference download is part of this static organization task.
