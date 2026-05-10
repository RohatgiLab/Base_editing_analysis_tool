### Base-editing screen analysis pipeline
A modular R-based pipeline for processing CRISPR gRNA base-editing screening data from raw FASTQ files through alignment, 
read counting, and statistical analysis (z-scores, p-values, and FDR), using control-guide log fold-changes as the null distribution.

This repository is designed to support reproducible, end-to-end analysis of pooled CRISPR base-editing screens, 
enabling standardized processing from sequencing reads to statistically inferred guide-level effects.

**An accompanying `.Rmd` file is provided to demonstrate a full example workflow, allowing users to run the pipeline step-by-step on example data and reproduce the analysis from raw sequencing reads to final statistical outputs.**
---

## Overview

This pipeline performs:

### Data Processing : fastq_to_readcounts_pipeline.R
- FASTQ preprocessing (trimming, motif-based filtering)
- gRNA reference index construction
- Alignment to gRNA library
- BAM sorting and indexing
- Read counting per gRNA
- Merging sample-level count matrices

### Statistical Analysis : readcounts_to_pvalues.R
- RPM normalization
- Replicate averaging
- Log2 transformation
- Log fold-change calculation
- Outlier filtering (IQR-based)
- Z-score normalization (control-based)
- p-value estimation
- Benjamini–Hochberg FDR correction

---


---

## Installation

Install required R packages:

```r
install.packages(c(
  "tidyverse",
  "readxl",
  "openxlsx",
  "stringr",
  "data.table"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Rsubread",
  "Biostrings",
  "GenomicAlignments",
  "GenomicFeatures",
  "QuasR",
  "Rsamtools",
  "ShortRead"
))


