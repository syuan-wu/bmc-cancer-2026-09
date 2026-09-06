# Multi-parameter Integration and Iterative Gene Selection for a ccRCC Metastasis-associated Signature

This repository contains the source code implementing the core computational workflow used in the study:

**“Multi-parameter integration algorithm via iterative sampling for identifying a metastasis-associated gene signature in clear cell renal cell carcinoma”**

The workflow integrates multi-parameter gene-level information, hierarchical clustering, iterative gene-combination sampling, logistic-regression-based feature selection, and gene-signature construction.

The code is organized into three independent R packages:

1. **MPICluster** – multi-parameter preprocessing, hierarchical clustering, and silhouette analysis
2. **IterGeneSelect** – iterative random gene-combination sampling and selection-frequency analysis
3. **SignatureModel** – logistic-regression-based feature selection, gene-combination search, final signature construction, and single-gene deletion analysis

The packages were designed to be **dataset-independent** so that the computational procedures can be inspected and reused with appropriately formatted datasets.

---

## Repository Structure

```text
.
├── MPICluster/
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   └── R/
│
├── IterGeneSelect/
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   └── R/
│
├── SignatureModel/
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   └── R/
│
├── full_pipeline_example.R
├── REPRODUCIBILITY_NOTES.md
└── README.md
