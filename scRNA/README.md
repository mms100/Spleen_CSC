# <img src="https://img.icons8.com/color/48/spleen.png" width="24" height="24" alt="spleen"/> Spleen CSC Single-Cell Analysis Toolkit
 This project is split into modular Nextflow pipelines designed to handle data processing from differential expression through downstream biological pathway validation.

---

##  Global Requirements

* **Nextflow Engine:** Both workflows require **Nextflow version `24.04.2.5914` or later** to run successfully.
* **Reproducibility:** Individual R script components and custom automation assets are housed inside the `./scripts` directory of each respective pipeline.

---

##  Project Architecture & Workflow

To run the complete analysis, follow the steps sequentially across the sub-modules:

### 1.  [DEG Pipeline](./scRNA-seq_nextflow_DEG_pipeline_nf)
The initial phase responsible for quality control, filtering, processing via Seurat, and identifying Differentially Expressed Genes (DEGs).
* *Go to this subdirectory to see raw data setup and execution scripts.*

### 2.  [Pathway Enrichment Pipeline](./Pathway_enrichment_pipeline_nf)
The downstream phase that ingests the metrics from the DEG step to perform pathway expansion and functional analysis.
* *Requires configuring the Zenodo mouse database reference in line 91.*


### 3.  [Cell-Cell communication]
Using Liana and crosstalkR to show the differential interactivity with bone marrow niche during fibrosis

