## TNBC Multi-Omics Target Discovery Pipeline Overview
This diagram summarizes the computational pipeline used to identify therapeutic targets in TNBC using TCGA transcriptomic data, network biology, and machine learning.
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
