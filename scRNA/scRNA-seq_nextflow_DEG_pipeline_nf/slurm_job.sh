#!/bin/bash
#SBATCH --job-name=DEG_pipeline     # Job name
#SBATCH --mail-user=email@email.com    # Where to send mail
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=path/to/scRNA-seq_nextflow_DEG_pipeline/log/output.%J.%x.txt
### Time to execute
#SBATCH --time=96:00:00

### amount of memory
#SBATCH --mem-per-cpu=8G

### amount of cores
#SBATCH --cpus-per-task=2

# Load necessary modules if you have on your HPC
# this modules have the following packages installed (MAST: 1.16.0, SingleCellExperiment: 1.28.1, Seurat: 4.4.0, SeuratObject: 4.1.4, dplyr: 1.2.1, optparse: 1.8.2, ggplot2: 4.0.3, tidyverse: 2.0.0, scales: 1.4.0,EnhancedVolcano: 1.24.0)
module load scRNA/1.0.4
module load R/4.4.1


# Run pipeline without batch
//path/to/executor/nextflow \
    run path/to/scRNA-seq_nextflow_DEG_pipeline/main.nf \
    --results_dir "path/to/scRNA-seq_nextflow_DEG_pipeline/output_WO_batch/"  \
    --object "path/to/scRNAseq_mouse_Thpo_EV.Rds" \
    --cond1 "ThPO" \
    --cond2 "EV" \
    --annotation "ECs_clustering" \
    --batch_colname "NULL"   \
    --cond_colname "stage" \
    --species "mouse"
    

    