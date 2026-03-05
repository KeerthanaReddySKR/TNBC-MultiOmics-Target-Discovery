# TNBC Multi-Omics Target Discovery – Pipeline Documentation

This document provides detailed documentation of the computational workflow used in the **TNBC Multi-Omics Target Discovery project**.  
The pipeline integrates transcriptomic analysis, mutation data, protein–protein interaction networks, and machine learning to identify candidate therapeutic targets in **Triple-Negative Breast Cancer (TNBC)** using TCGA data.

---

# Pipeline Overview

The computational pipeline follows the workflow below:

TCGA BRCA Dataset  
↓  
Expression Import & Audit  
↓  
Sample Quality Control  
↓  
TNBC Cohort Selection  
↓  
Differential Expression Analysis (limma)  
↓  
Gene Set Enrichment Analysis (MSigDB Hallmark)  
↓  
Mutation Data Integration  
↓  
STRING Protein–Protein Interaction Network  
↓  
Network Sensitivity Analysis  
↓  
Machine Learning Target Ranking  
↓  
STRING Gene Cutoff Determination  
↓  
Final Therapeutic Target Identification  

---

# Dataset

**Source:** The Cancer Genome Atlas (TCGA)

Datasets used in this project:

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

Operations performed:

- Load expression matrix from TCGA Xena
- Convert expression table into a gene × sample matrix
- Verify TCGA barcode structure
- Detect duplicate genes or samples
- Load clinical metadata
- Identify overlapping samples between expression and clinical data
- Subset expression matrix to valid samples

Output:

data_processed/expr_clean.rds

This file contains the cleaned expression matrix used for downstream analysis.

---

## 2. TNBC Label Construction

Script:

scripts/02_define_TNBC_and_build_labels.R

TNBC samples are defined using receptor status information from clinical annotations.

Criteria used:

- ER negative
- PR negative
- HER2 negative

Processing steps include:

- Extract receptor status from clinical table
- Normalize receptor status labels
- Identify TNBC samples
- Remove samples with missing receptor information
- Generate TNBC classification labels

Outputs:

data_processed/BRCA_TNBC_labels.rds  
data_processed/BRCA_TNBC_labels.tsv  

---

## 3. Sample Quality Control Summary

Script:

scripts/00_QC_sample_flow.R

This script summarizes how raw TCGA samples were filtered during preprocessing.

The script:

- Loads the raw expression dataset
- Loads the cleaned dataset
- Compares number of genes and samples
- Reports how many samples were removed during QC

Output:

results/QC_sample_flow_summary.tsv

This file summarizes:

- Total downloaded samples
- Final samples used in analysis
- Samples removed during quality control

---

## 4. Differential Gene Expression Analysis

Script:

scripts/05_differential_expression_limma.R

Differential expression analysis identifies genes significantly dysregulated in TNBC.

Method used:

limma

Outputs include:

- Log fold change
- Adjusted p-values
- Ranked gene lists

Results stored in:

results/DEG/

---

## 5. Gene Set Enrichment Analysis

Script:

scripts/06_GSEA_hallmark_analysis.R

Gene Set Enrichment Analysis (GSEA) identifies biological pathways associated with TNBC gene expression.

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
- Mutated in TNBC

Combining expression and mutation information improves biological relevance for downstream analysis.

---

## 7. STRING Protein Interaction Network

Script:

scripts/08_STRING_network_analysis.R

A protein–protein interaction network is constructed using the STRING database.

Purpose:

- Identify functional interactions among dysregulated genes
- Detect network modules related to TNBC biology

Results stored in:

results/STRING/

---

## 8. Network Sensitivity Analysis

Script:

scripts/08b_network_sensitivity_analysis.R

Network sensitivity analysis evaluates the robustness of the STRING interaction network.

This step helps determine:

- Network stability
- Key regulatory genes
- Highly connected nodes in the interaction network

---

## 9. Machine Learning Target Ranking

Script:

scripts/10_network_guided_ML_target_ranking.R

Machine learning models are used to prioritize candidate genes.

Features used include:

- Differential expression magnitude
- Pathway enrichment signals
- Mutation status
- Network connectivity

Model performance is evaluated using **nested cross-validation**.

Results stored in:

results/ML/

---

## 10. STRING Gene Cutoff Determination

Script:

scripts/11_define_STRING_gene_cutoff.R

To determine an objective cutoff for STRING network analysis, cumulative contribution of normalized gene scores is calculated.

Steps include:

- Normalize gene importance scores
- Compute cumulative contribution
- Identify gene thresholds representing **85% and 90% of signal**

Outputs generated:

Genes_for_STRING_85pct.tsv  
Genes_for_STRING_90pct.tsv  
Cumulative_score_curve.png  

---

# Final Outputs

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

environment/software_environment.md

3. Execute scripts sequentially from the **scripts/** directory.

---

# Author

Y. Showri Keerthana  
B.Tech Bioinformatics  
VFSTR University, India

Research Interests:

- Computational Biology
- Cancer Genomics
- Network Biology
- Machine Learning for Biomedical Discovery
