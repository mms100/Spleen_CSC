#!/bin/bash
#SBATCH --job-name=pathway_Thpo_EV     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=/path/to/logs/output.%J.%x.txt
#SBATCH --mail-user=email@email.com   # Where to send mail	
### Time to execute
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=30G

### OpenMP threads
#SBATCH --cpus-per-task=10



#running the pathway analysis

#load needed environment
# this modules have the following packages installed (MAST: 1.16.0, SingleCellExperiment: 1.28.1, Seurat: 4.4.0, SeuratObject: 4.1.4, dplyr: 1.2.1, optparse: 1.8.2, ggplot2: 4.0.3, tidyverse: 2.0.0, scales: 1.4.0,EnhancedVolcano: 1.24.0)
# to recreate the environment run this: conda env create -f environment.yml
source //path/to/miniconda/bin/activate
conda activate Decopler
module load scRNA/1.0.4
module load R/4.4.1 

//path/to/nextflow \
    run pathway_main.nf \
    --inputdir "/path/to/DEGs_tables/" \
    --cond1 "ThPO" \
    --cond2 "EV" \
    --species "mouse" 
