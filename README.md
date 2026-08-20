# cfRNA_microbiome
Integrative analysis of plasma cell-free RNA and gut microbiome reveal novel molecular signatures in GDM

# Multi-Omics Visualization Pipeline for GDM Research

This repository provides a set of self contained R scripts to reproduce key figures from a multi omics study of gestational diabetes mellitus (GDM). The pipeline covers cell free RNA (cfRNA), gut microbiome (KO), and metabolomics data, including differential expression, pathway enrichment, WGCNA, MOFA integration, and machine learning prediction.

All scripts use simulated data to demonstrate the plotting logic. They are ready to run without any external input files and can be easily adapted to real datasets.

---

## Repository Contents

| File | Description |
|------|-------------|
| Data processing.r | Simulated preprocessing: low expression filtering, variance selection, log2 transformation, ComBat batch correction, MOFA input construction, and sample anonymization. |
| fig2.r | Figure 2   Volcano plot (DEGs) + GSEA bubble plot (KEGG pathways across trimesters). |
| fig3.r | Figure 3   Donut charts (cell type origins), trophoblast boxplots, WGCNA module enrichment heatmap, and module trait correlation heatmap. |
| fig4.r | Figure 4   Species dotplot, circular GSEA fan plot for microbiome KOs, and pathway specific boxplots. |
| fig5.r | Figure 5   MOFA variance explained barplot, factor trait correlation heatmap, and back to back lollipop plot for miRNA target enrichment. |
| fig6.r | Figure 6   ROC curves for GDM prediction models (cfRNA, cfRNA+BMI, cfRNA+KO) and BMI z score prediction scatter/selection plots. |
| Metabolomics analysis.r | Metabolite differential testing (Wilcoxon) for BCAA, SCFA, and bile acids, plus volcano and boxplot visualizations. |
| Sensitivity analysis.r | Sensitivity metrics (Jaccard, Spearman, Kappa) and scatter plots comparing alternative covariate models. |
| README.md | This file. |

---

## Software and Algorithms

| Software / Package | Version | Reference / URL |
|---------------------|---------|-----------------|
| R | 4.2.2 | [https://www.r-project.org/](https://www.r-project.org/) |
| limma | 3.54.2 | [https://bioinf.wehi.edu.au/limma/](https://bioinf.wehi.edu.au/limma/) |
| ggplot2 | 4.0.3 | [https://ggplot2.tidyverse.org/](https://ggplot2.tidyverse.org/) |
| ggraph | 2.1.0 | [https://github.com/thomasp85/ggraph](https://github.com/thomasp85/ggraph) |
| lme4 | 1.1 31 | [https://github.com/lme4/lme4](https://github.com/lme4/lme4) |
| lmerTest | 3.1 3 | [https://github.com/runehaubo/lmerTestR](https://github.com/runehaubo/lmerTestR) |
| clusterProfiler | 4.6.0 | Yu, G. et al.[64,65]   [https://github.com/YuLab-SMU/clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) |
| MOFA2 | (latest) | Argelaguet, R. et al.[66]   [https://github.com/bioFAM/MOFA2](https://github.com/bioFAM/MOFA2) |
| scikit learn (sklearn) | 1.1.2 | Pedregosa, F. et al.[67]   [https://scikit-learn.org/](https://scikit-learn.org/) |
| TensorFlow | 2.9.1 | Martin, A.[68]   [https://www.tensorflow.org/](https://www.tensorflow.org/) |

> Note: The scripts in this repository primarily use R packages. Python packages (scikit learn, TensorFlow) are listed for completeness as they were used in the original study, but the provided R scripts rely only on R dependencies (the figure generation does not require Python). The R packages can be installed as follows:

r
# CRAN packages
install.packages(c(
  "ggplot2", "ggrepel", "dplyr", "tidyr", "patchwork",
  "MatrixGenerics", "pROC", "lme4", "lmerTest"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "limma", "ComplexHeatmap", "circlize", "sva", "clusterProfiler"
))


---

## Usage

### Run a single script
bash
Rscript fig2.r


### Run all figure scripts
bash
for f in fig.r; do Rscript $f; done


### Adapt to your own data
1. Replace the mock data generation blocks (e.g., set.seed() + rnorm()/runif()) with your actual data frames.
2. Ensure column names match those expected in the plotting code (e.g., logFC, p.adjust, NES, pathway).
3. Adjust thresholds (e.g., logFC > 0.15, p.adjust < 0.05) to fit your analysis.

---

## Output

| Script | Output |
|--------|--------|
| fig2.r   fig6.r | Displays plots in the active graphics device. |
| fig4.r | Also saves Figure4B_fan.pdf (circular GSEA plot) in the working directory. |
| Metabolomics analysis.r | Displays volcano and boxplots; writes metabolite_results.csv (mock results). |
| Sensitivity analysis.r | Displays scatter and bar plots; prints metric table to console. |

To save any plot as a PDF, uncomment the ggsave() lines at the end of the scripts.

---

## Notes

- All scripts are fully self contained and require no external data files.
- The circular GSEA plot in fig4.r uses circlize and is saved as a PDF because it cannot be embedded directly in a ggplot layout.
- The scripts are designed for reproducibility and educational purposes; they illustrate the visualisation pipeline without exposing real patient data.

---

## License

This project is licensed under the MIT License   see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use these scripts in your work, please cite our manuscript (citation details will be added upon publication).

---

## Contact

For questions or suggestions, please open an issue on this repository or contact the corresponding author.
