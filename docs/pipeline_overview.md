# TNBC Multi-Omics Target Discovery – Pipeline Documentation

This document provides detailed documentation of the computational workflow used in the **TNBC Multi-Omics Target Discovery project**.  
The pipeline integrates transcriptomic analysis, mutation data, protein–protein interaction networks, and machine learning to identify candidate therapeutic targets in **Triple-Negative Breast Cancer (TNBC)** using TCGA data.

---

# Pipeline Overview

The computational pipeline consists of the following stages:

1. TCGA BRCA Dataset Acquisition  
2. Expression Data Import and Audit  
3. Sample Quality Control  
4. TNBC Cohort Selection  
5. Differential Expression Analysis (limma)  
6. Gene Set Enrichment Analysis (MSigDB Hallmark)  
7. Mutation Data Integration  
8. STRING Protein–Protein Interaction Network Construction  
9. Network Sensitivity Analysis  
10. Machine Learning Target Ranking  
11. STRING Gene Cutoff Determination  
12. Final Therapeutic Target Identification  

---

# Dataset

**Source:** The Cancer Genome Atlas (TCGA)

Datasets used in this project include:

- TCGA-BRCA RNA-seq transcriptomic data  
- Clinical metadata  
- Somatic mutation data  

The dataset was filtered to extract **Triple-Negative Breast Cancer (TNBC)** samples based on receptor status.

---

# Pipeline Steps

## 1. Expression Data Import and Audit

Script:

scripts/01_import_and_audit_expression.R

This step imports the TCGA RNA-seq expression matrix and verifies data integrity.

Key operations:

- Load expression matrix from TCGA Xena
- Convert expression table into gene × sample matrix
- Verify TCGA barcode format
- Detect duplicate genes or samples
- Load clinical metadata
- Identify overlapping samples between expression and clinical datasets
- Subset expression matrix to valid samples

Output generated:

data_processed/expr_clean.rds

This file contains the cleaned expression matrix used for downstream analysis.

---

## 2. TNBC Label Construction

Script:

scripts/02_define_TNBC_and_build_labels.R

TNBC samples are defined based on receptor status from clinical annotations.

Criteria used:

- ER negative  
- PR negative  
- HER2 negative  

Processing steps include:

- Extract receptor status from clinical data
- Normalize receptor status labels
- Identify TNBC samples
- Remove samples with missing receptor information
- Construct TNBC classification labels

Outputs generated:

data_processed/BRCA_TNBC_labels.rds  
data_processed/BRCA_TNBC_labels.tsv  

---

## 3. Sample Quality Control Summary

Script:

scripts/00_QC_sample_flow.R

This script summarizes how raw TCGA samples were filtered during preprocessing.

The script performs:

- Loading of raw expression dataset
- Loading of cleaned dataset
- Comparison of sample counts
- Reporting number of removed samples during QC

Output generated:

results/QC_sample_flow_summary.tsv

This file summarizes:

- Downloaded samples
- Final samples used
- Samples removed during quality control

---

## 4. Differential Gene Expression Analysis

Script:

scripts/05_differential_expression_limma.R

Differential expression analysis identifies genes significantly dysregulated in TNBC.

Method used:

limma

Outputs include:

- Log fold change values
- Adjusted p-values
- Ranked gene lists

Results stored in:

results/DEG/

---

## 5. Gene Set Enrichment Analysis

Script:

scripts/06_GSEA_hallmark_analysis.R

Gene Set Enrichment Analysis (GSEA) identifies biological pathways associated with TNBC gene expression patterns.

Database used:

MSigDB Hallmark pathways

Outputs include:

- Enriched biological pathways
- Normalized enrichment scores
- Statistical significance values

Results stored in:

results/GSEA/

---

## 6. Mutation Data Integration

Script:

scripts/07_mutation_integration.R

Somatic mutation data from TCGA are integrated with transcriptomic results.

This step helps identify genes that are both:

- Differentially expressed
- Mutationally altered in TNBC

Combining expression and mutation information improves biological interpretation.

---

## 7. STRING Protein Interaction Network

Script:

scripts/08_STRING_network_analysis.R

Protein–protein interaction networks are constructed using the **STRING database**.

Purpose of this step:

- Identify functional interactions among dysregulated genes
- Detect network modules associated with TNBC biology

Results stored in:

results/STRING/

---

## 8. Network Sensitivity Analysis

Script:

scripts/08b_network_sensitivity_analysis.R

Network sensitivity analysis evaluates the robustness of the protein interaction network.

This step identifies:

- Highly connected genes
- Influential nodes
- Potential regulatory hubs

---

## 9. Machine Learning Target Ranking

Script:

scripts/10_network_guided_ML_target_ranking.R

Machine learning models are used to prioritize candidate therapeutic targets.

Features used for ranking include:

- Differential expression magnitude
- Pathway enrichment information
- Mutation status
- Network connectivity

Models are evaluated using **nested cross-validation**.

Results stored in:

results/ML/

---

## 10. STRING Gene Cutoff Determination

Script:

scripts/11_define_STRING_gene_cutoff.R

An objective cutoff for STRING network analysis is determined using cumulative contribution of normalized gene scores.

Steps include:

- Normalize gene importance scores
- Compute cumulative contribution
- Identify thresholds covering **85% and 90% of signal**

Outputs generated:

Genes_for_STRING_85pct.tsv  
Genes_for_STRING_90pct.tsv  
Cumulative_score_curve.png  

---

# Pipeline Outputs

The pipeline generates several intermediate and final outputs including:

- Quality control reports
- Differential expression results
- Pathway enrichment results
- Mutation comparison tables
- Protein interaction networks
- Machine learning ranking scores
- Candidate therapeutic targets

All outputs are organized within the **results/** directory.

---

# Reproducibility

To reproduce the analysis:

1. Download TCGA BRCA RNA-seq data.  
2. Install required R packages listed in:
   **environment/software_environment.md**
3. Execute scripts sequentially from the **scripts/** directory.
