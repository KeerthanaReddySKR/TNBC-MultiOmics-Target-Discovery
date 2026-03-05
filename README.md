# TNBC Multi-Omics Target Discovery

A network-guided computational pipeline for identifying potential therapeutic targets in **Triple-Negative Breast Cancer (TNBC)** by integrating **TCGA transcriptomic data**, **protein–protein interaction networks**, and **machine learning-based prioritization**.

The pipeline combines transcriptomic analysis, pathway enrichment, mutation integration, network biology, and predictive modeling to systematically prioritize candidate therapeutic targets.

---
## Table of Contents

* [Pipeline Overview](#pipeline-overview)
* [Dataset](#dataset)
* [Computational Workflow](#computational-workflow)
* [Repository Structure](#repository-structure)
* [Key Methods and Tools](#key-methods-and-tools)
* [Results](#results)
* [Reproducibility](#reproducibility)
* [Author](#author)

# Pipeline Overview

This pipeline processes TCGA breast cancer transcriptomic data to identify TNBC-specific gene expression patterns and prioritize therapeutic targets through network-guided machine learning.

## Computational Workflow

```
TCGA BRCA Dataset
        ↓
Sample QC & Filtering
        ↓
TNBC Cohort Selection
        ↓
Differential Expression (limma)
        ↓
GSEA Pathway Analysis
        ↓
Mutation Integration
        ↓
STRING PPI Network
        ↓
Network Sensitivity Analysis
        ↓
Machine Learning Target Ranking
        ↓
Final Therapeutic Targets
```

For a step-by-step explanation of each stage of the analysis pipeline, see the detailed documentation:

[Pipeline Step Documentation](docs/pipeline_overview.md)

---

# Dataset

**Source:** The Cancer Genome Atlas (TCGA)

Datasets used in this project:

* **TCGA-BRCA RNA-seq transcriptomic data**
* **Clinical metadata**
* **Somatic mutation data**

The dataset was filtered to extract **TNBC patient samples** based on receptor status.

---

# Computational Workflow

The analysis pipeline consists of the following stages:

1. Sample quality control and filtering
2. TNBC cohort selection
3. Differential gene expression analysis using **limma**
4. Gene set enrichment analysis using **MSigDB Hallmark pathways**
5. Mutation data integration
6. Construction of **STRING protein–protein interaction networks**
7. Network sensitivity analysis
8. Machine learning-based target ranking
9. Final candidate therapeutic target identification

---

# Repository Structure

```
TNBC-MultiOmics-Target-Discovery
│
├── docs
│   └── pipeline_overview.md
│
├── environment
│   └── software_environment.md
│
├── scripts
│   ├── 00_QC_sample_flow.R
│   ├── 01_import_and_audit_expression.R
│   ├── 02_define_TNBC_labels.R
│   ├── 03_build_final_cohorts.R
│   ├── 04_transcriptome_QC_PCA.R
│   ├── 05_differential_expression_limma.R
│   ├── 06_GSEA_hallmark_analysis.R
│   ├── 07_mutation_integration.R
│   ├── 08_STRING_network_analysis.R
│   ├── 08b_network_sensitivity_analysis.R
│   ├── 09_define_STRING_gene_sets.R
│   ├── 10_network_guided_ML_target_ranking.R
│   ├── 11_define_STRING_gene_cutoff.R
│   ├── 12_nested_cross_validation.R
│   ├── 13_final_target_panel.R
│   └── 14_final_target_intersection.R
│
├── results
│   ├── QC
│   ├── DEG
│   ├── GSEA
│   ├── STRING
│   ├── ML
│   └── targets
│
├── figures
│   └── pipeline_overview.md
│
├── README.md
└── .gitignore
```

---

# Key Methods and Tools

The pipeline integrates multiple computational biology tools and resources:

* **R programming language**
* **limma** for differential expression analysis
* **MSigDB Hallmark pathways** for gene set enrichment analysis
* **STRING database** for protein–protein interaction networks
* **Machine learning models** for therapeutic target prioritization

---

# Results

The pipeline identifies candidate genes that may serve as **therapeutic targets in Triple-Negative Breast Cancer (TNBC)** by integrating transcriptomic signals with protein interaction networks and predictive modeling.

Intermediate outputs and final results are organized within the **results/** directory.
## Results

### TCGA Data Acquisition

RNA-seq gene expression data for the TCGA Breast Cancer (BRCA) cohort was downloaded using the TCGAbiolinks package.

Dataset characteristics:

- Total samples: 1231
- Data type: Gene Expression Quantification
- Workflow: STAR – Counts
- Reference genome: hg38

### Expression Matrix Construction

Raw STAR count files were processed and merged into a unified expression matrix.

Final matrix dimensions:

- Genes: 60,664
- Samples: 1,231

Each column corresponds to a TCGA sample, and each row represents a gene feature.

### Data Structure

The expression matrix was converted into a Bioconductor `SummarizedExperiment` object for compatibility with downstream genomic analysis tools.

Object summary:

Class: SummarizedExperiment  
Dimensions: 60,664 genes × 1,231 samples  
Assay: raw STAR count matrix

### Final Dataset

The processed dataset was saved in R serialized format:

data_raw/TCGA_BRCA_STAR_counts_FINAL.rds

This file serves as the primary input for downstream analyses including TNBC cohort selection, differential expression analysis, pathway enrichment, network analysis, and machine learning–based target prioritization.

---

# Reproducibility

To reproduce the analysis:

1. Download TCGA BRCA RNA-seq data.
2. Install the required R packages listed in:

```
environment/software_environment.md
```

3. Execute the scripts sequentially from the **scripts/** directory.

---

# Author

**Y. Showri Keerthana**
B.Tech Bioinformatics
VFSTR University, India

### Research Interests

* Computational Biology
* Cancer Genomics
* Network Biology
* Machine Learning for Biomedical Discovery
