#  Pathway Enrichment Analysis Pipeline

This repository contains the Nextflow pipeline for downstream **Pathway Enrichment Analysis**. 

>  **Prerequisite:** This pipeline is dependent on, and runs using, the output results generated from the [DEG Analysis Pipeline](../scRNA-seq_nextflow_DEG_pipeline_nf).

---

##  Getting Started

### 1. Environment Setup
To replicate the exact environment and install all the necessary R/Bioconductor dependencies required for this analysis, run the following command using Conda:

```bash
conda env create -f environment.yml
