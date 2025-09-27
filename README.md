# RNA-seq Analysis for "Targeting KRAS inhibitor-resistant pancreatic cancer with a M1C antibody-drug conjugate"
This repository contains the R code used to generate the figures for our manuscript: Hiroki Ozawa. Kazuki Takahashi.. Tomoki Motegi. Justine Jacobi. Atrayee Bhattarchya. Keisuke Shigeta. Shinkichi Takamori. Mai Onishi. Zhou Luan. Kazumasa Fukuda. Akihisa Ueno. Minoru Kitago. Hirofumi Kawakubo. Yuko Kitagawa. Dove Keith. John Muschler. Rosalie Sears. Mark D. Long, Joseph D. Mancias, Andrew J. Aguirre and Donald Kufe*: "Targeting KRAS inhibitor-resistant pancreatic cancer with a M1C antibody-drug conjugate"  
It covers preprocessing of count data, differential expression, batch correction, consensus clustering, gene-set analyses (Hallmark, stemness, basal/classical, mucinous differentiation), and figure generation (volcano, heatmaps, scatter/violin).

> **Data note**  
> Raw FASTQ were processed with **STAR–RSEM** (external to this repo). Downstream analyses here start from gene-level count/CPM tables.


---

## Requirements

- **R** ≥ 4.2
- CRAN/Bioconductor packages:

```r
# Bioconductor core
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "TCC", "sva", "genefilter", "clusterProfiler",
  "enrichplot", "ComplexHeatmap"
))

# CRAN
install.packages(c(
  "tidyverse", "ggplot2", "ggrepel", "ConsensusClusterPlus",
  "circlize", "grid", "msigdbr", "multcomp", "ggpp", "ggpmisc"
))

