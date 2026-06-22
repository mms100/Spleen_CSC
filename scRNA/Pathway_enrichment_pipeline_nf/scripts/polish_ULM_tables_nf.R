# Load required libraries
library(readr)
library(dplyr)
library(writexl)
library(readxl)

# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) return(NA)
  args[idx + 1]
}

ulm_tables_dir <- get_arg('--ulm_tables_dir')
gene_annotation_path <- get_arg('--gene_annotation')
cond1 <- get_arg('--cond1')
cond2 <- get_arg('--cond2')
csv_path_base <- get_arg('--csv_path_base')
if (is.na(csv_path_base)) {
  csv_path_base <- file.path(dirname(csv_path_base), 'tables')
}

if (is.na(ulm_tables_dir) || is.na(gene_annotation_path) || is.na(cond1) || is.na(cond2)) {
  stop('Usage: Rscript polish_ULM_tables_nf.R --ulm_tables_dir <ULM_tables_dir> --gene_annotation <gene_annotation.csv> --cond1 <cond1> --cond2 <cond2> [--csv_path_base <csv_path>]')
}

comparison <- paste0(cond1, 'vs', cond2)

collections_to_process <- c('GO', 'Hallmark', 'KEGG', 'Reactome')

gene_annotation <- read.csv(gene_annotation_path, row.names = 1)
gene_annotation <- gene_annotation[, c('geneset', 'genesymbol', 'collection')]

for (collection in collections_to_process) {
  cat('Processing collection:', collection, '\n')
  acts_path <- file.path(ulm_tables_dir, collection, 'acts.csv')
  pval_path <- file.path(ulm_tables_dir, collection, 'pvals.csv')
  if (!file.exists(acts_path) | !file.exists(pval_path)) {
    cat('  Skipping', collection, '- acts.csv or pvals.csv not found.\n')
    next
  }
  acts_df <- read.csv(acts_path, row.names = 1)
  pval_df <- read.csv(pval_path, row.names = 1)
  sheets_list_filtered <- list()
  sheets_list_unfiltered <- list()
  for (col in colnames(acts_df)) {
    if (col != 'terms') {
      combined_df <- data.frame(
        geneset = rownames(acts_df),
        ulm_score = acts_df[[col]],
        pval = as.numeric(pval_df[[col]])
      )
      filtered_df <- combined_df %>% filter(pval <= 0.05)
      sheets_list_filtered[[col]] <- filtered_df
      sheets_list_unfiltered[[col]] <- combined_df
    }
  }
  excel_file_filtered <- paste0("combined_", collection, "_filtered.xlsx")
  excel_file_unfiltered <- paste0("combined_", collection, "_unfiltered.xlsx")
  write_xlsx(sheets_list_filtered, excel_file_filtered)
  cat('  Wrote filtered data to', excel_file_filtered, '\n')
  write_xlsx(sheets_list_unfiltered, excel_file_unfiltered)
  cat('  Wrote unfiltered data to', excel_file_unfiltered, '\n')

  # Subset gene annotation
  if (tolower(collection) == 'go') {
    gene_annotation_sub <- gene_annotation[grepl('^GOBP_|^GOCC_|^GOMF_', gene_annotation$geneset), ]
    gene_annotation_sub$geneset <- gsub('^GOBP_|^GOCC_|^GOMF_', '', gene_annotation_sub$geneset)
  } else if (tolower(collection) == 'hallmark') {
    gene_annotation_sub <- gene_annotation[grepl('^HALLMARK_', gene_annotation$geneset), ]
    gene_annotation_sub$geneset <- gsub('^HALLMARK_', '', gene_annotation_sub$geneset)
  } else if (tolower(collection) == 'kegg') {
    gene_annotation_sub <- gene_annotation[grepl('^KEGG_', gene_annotation$geneset), ]
    gene_annotation_sub$geneset <- gsub('^KEGG_', '', gene_annotation_sub$geneset)
  } else if (tolower(collection) == 'reactome') {
    gene_annotation_sub <- gene_annotation[grepl('^REACTOME_', gene_annotation$geneset), ]
    gene_annotation_sub$geneset <- gsub('^REACTOME_', '', gene_annotation_sub$geneset)
  } else {
    gene_annotation_sub <- gene_annotation[0, ]
  }

  csv_files <- list.files(path = csv_path_base, pattern = '*.csv', full.names = TRUE)
  pattern <- paste0('^', comparison, '_|\\.csv$')
  csv_names <- gsub(pattern, '', basename(csv_files))

  sheet_names_filtered <- excel_sheets(excel_file_filtered)
  excel_data_filtered <- setNames(
    lapply(sheet_names_filtered, function(sheet) read_excel(excel_file_filtered, sheet = sheet)),
    sheet_names_filtered
  )
  updated_sheets_filtered <- list()
  for (i in seq_along(csv_files)) {
    csv_name <- csv_names[i]
    if (csv_name %in% sheet_names_filtered) {
      df <- read_csv(csv_files[i], show_col_types = FALSE)
      if (!'GeneID' %in% colnames(df)) {
        cat(sprintf('  Warning: No GeneID column in %s, skipping.\n', csv_files[i]))
        next
      }
      filtered_genes <- gene_annotation_sub %>% filter(genesymbol %in% df$GeneID)
      sheet_data <- excel_data_filtered[[csv_name]]
      colnames(sheet_data) <- c('geneset', 'ulm_score', 'pval')
      sheet_data$Matched_Genes <- sapply(sheet_data$geneset, function(gs) {
        matched_genes <- filtered_genes$genesymbol[filtered_genes$geneset == gs]
        if (length(matched_genes) == 0) return(NA)
        paste(matched_genes, collapse = ', ')
      })
      updated_sheets_filtered[[csv_name]] <- sheet_data
    }
  }
  write_xlsx(updated_sheets_filtered, excel_file_filtered)
  cat('  Annotated filtered Excel file for', collection, '\n')

  sheet_names_unfiltered <- excel_sheets(excel_file_unfiltered)
  excel_data_unfiltered <- setNames(
    lapply(sheet_names_unfiltered, function(sheet) read_excel(excel_file_unfiltered, sheet = sheet)),
    sheet_names_unfiltered
  )
  updated_sheets_unfiltered <- list()
  for (i in seq_along(csv_files)) {
    csv_name <- csv_names[i]
    if (csv_name %in% sheet_names_unfiltered) {
      df <- read_csv(csv_files[i], show_col_types = FALSE)
      if (!'GeneID' %in% colnames(df)) {
        cat(sprintf('  Warning: No GeneID column in %s, skipping.\n', csv_files[i]))
        next
      }
      filtered_genes <- gene_annotation_sub %>% filter(genesymbol %in% df$GeneID)
      sheet_data <- excel_data_unfiltered[[csv_name]]
      colnames(sheet_data) <- c('geneset', 'ulm_score', 'pval')
      sheet_data$Matched_Genes <- sapply(sheet_data$geneset, function(gs) {
        matched_genes <- filtered_genes$genesymbol[filtered_genes$geneset == gs]
        if (length(matched_genes) == 0) return(NA)
        paste(matched_genes, collapse = ', ')
      })
      updated_sheets_unfiltered[[csv_name]] <- sheet_data
    }
  }
  write_xlsx(updated_sheets_unfiltered, excel_file_unfiltered)
  cat('  Annotated unfiltered Excel file for', collection, '\n')
}
cat('Processing complete! All collections have been processed and saved.\n')


