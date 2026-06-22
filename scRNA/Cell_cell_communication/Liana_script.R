#load_packages

library(Seurat)
library(dplyr)
library(liana)
library(tibble)
library(archive)

load_object <- function(file_name){
 con <- archive::file_read(file = file_name)
 res <- readRDS(file = con)
 close(con)
 return(res)}

#load_object
spleen_object <- load_object("//path/to/scRNAseq_mouse_Thpo_EV.Rds")

#remove BCell as they only have 6 cells in TPO
light_spleen <- subset(spleen_object, subset = cell_type_final == "BCell" , invert = T)

#run the crosstalkR on the final annotation
outpath <- "//path/to/crosstalkR_analysis/"
pval_threshold <- 0.05

spleen_object_list <- SplitObject(light_spleen, split.by = "stage")

spleen_object_list$ThPO <- SetIdent(spleen_object_list$ThPO, value = "ECs_clustering")
spleen_object_list$EV <- SetIdent(spleen_object_list$EV, value = "ECs_clustering")

for (condition in names(spleen_object_list)) {
  sub_object <- spleen_object_list[[condition]]
  
  liana_results <- liana_wrap(sub_object, 
                              method = "cellphonedb",
                              resource = c("MouseConsensus"),
                              expr_prop = 0.1
  )
  
  liana_results <- liana_results %>%
    filter(pvalue < pval_threshold) %>%
    select(source, target, ligand, receptor.complex, lr.mean) %>%
    rename(gene_A = ligand, gene_B = receptor.complex, MeanLR = lr.mean)
  liana_results$type_gene_A = 'Ligand'
  liana_results$type_gene_B = 'Receptor'
  
  write.csv(liana_results, paste0(outpath, condition, "_LR.csv"))
}

