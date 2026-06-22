# Load necessary libraries
library(dendextend)
library(dplyr)
library(ComplexHeatmap)
library(optparse)
library(grid)

# parser for output and input directories
args <- commandArgs(trailingOnly = TRUE)
inputdir <- args[which(args == "--input_dir_6") + 1]
outdir <- args[which(args == "--outdir_2") + 1]

# Function to read CSV files from a given folder
read_csv_files_from_folder <- function(folder_path) {
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  subsubfolder_names <- basename(csv_files)
  dataframes <- lapply(csv_files, function(file) {
    df <- as.matrix(read.csv(file = file, check.names = FALSE, row.names = 1))
    return(df)
  })
  names(dataframes) <- subsubfolder_names
  return(dataframes)
}

# Function to read CSV files from all subfolders in the main folder
main_folder_to_list_of_lists <- function(main_folder_path) {
  subfolders <- list.dirs(path = main_folder_path, full.names = TRUE, recursive = FALSE)
  subfolder_names <- basename(subfolders)
  all_dataframes <- lapply(subfolders, read_csv_files_from_folder)
  names(all_dataframes) <- subfolder_names
  return(all_dataframes)
}

main_folder_path <- inputdir
all_csv_data_1 <- main_folder_to_list_of_lists(main_folder_path)

process_and_plot_heatmap <- function(data_list, heatmap_title) {
  # extract a matrix of the combined enrichment score for the heatmap
  combined_df_heatmap <- as.matrix(data_list$acts)
  combined_df_heatmap[is.na(combined_df_heatmap)] <- 0
  combined_df_heatmap[is.nan(combined_df_heatmap)] <- 0
  combined_df_heatmap[is.infinite(combined_df_heatmap)] <- 0
  # Scale the matrix for visualization
  combined_df_heatmap <- scale(combined_df_heatmap)
  
  # Make a matrix for the significance tags
  tags_df_heatmap <- as.matrix(data_list$pvals)
  tags_df_heatmap[is.na(tags_df_heatmap)] <- 1
  tags_df_heatmap[is.nan(tags_df_heatmap)] <- 0
  tags_df_heatmap[is.infinite(tags_df_heatmap)] <- 0

  # If more than 50 rows, keep only those containing the top 50 values in the matrix
  #if (nrow(combined_df_heatmap) > 50) {
    #mat_df <- as.data.frame(as.table(combined_df_heatmap))
    #colnames(mat_df) <- c("Row", "Column", "Value")
    #top_50 <- mat_df[order(-abs(mat_df$Value)), ][1:100, ]
    #top_50_rows <- unique(top_50$Row)
    #combined_df_heatmap <- combined_df_heatmap[rownames(combined_df_heatmap) %in% top_50_rows, , drop=FALSE]
    #tags_df_heatmap <- tags_df_heatmap[rownames(tags_df_heatmap) %in% top_50_rows, , drop=FALSE]
  #}



  if (nrow(combined_df_heatmap) > 50) {
  mat_df <- as.data.frame(as.table(combined_df_heatmap))
  colnames(mat_df) <- c("Row", "Column", "Value")
  # Order by absolute value, descending
  mat_df <- mat_df[order(-abs(mat_df$Value)), ]

  # Collect unique rows in order of top values until we have 50
  unique_rows <- character(0)
  for (row in mat_df$Row) {
    if (!(row %in% unique_rows)) {
      unique_rows <- c(unique_rows, row)
      if (length(unique_rows) == 50) break
    }
  }

  combined_df_heatmap <- combined_df_heatmap[rownames(combined_df_heatmap) %in% unique_rows, , drop=FALSE]
  tags_df_heatmap <- tags_df_heatmap[rownames(tags_df_heatmap) %in% unique_rows, , drop=FALSE]
  }
  #start the function for the dendrogram 
  customize_dendrogram <- function(hc, lwd = 2) {
    dend <- as.dendrogram(hc)
    dend <- dendextend::set(dend, "branches_lwd", lwd)
    return(dend)
  }
  
  # Generate the heatmap with customized dendrogram
  hc_row <- hclust(dist(combined_df_heatmap))
  hc_col <- hclust(dist(t(combined_df_heatmap)))
  row_dend <- customize_dendrogram(hc_row, lwd = 0.5)
  col_dend <- customize_dendrogram(hc_col, lwd = 0.5)
  
  heatmap <- Heatmap(combined_df_heatmap,
                     width = ncol(combined_df_heatmap)*unit(2, "mm"), 
                     height = nrow(combined_df_heatmap)*unit(1, "mm"),
                     heatmap_legend_param = list(
                       legend_direction = "horizontal",
                       title_gp = gpar(fontsize = 5, fontface = "bold"),
                       labels_gp = gpar(fontsize = 2)
                     ),
		     cell_fun = function(j, i, x, y, w, h, fill) {
			     pval <- tags_df_heatmap[i, j]
			     if (pval < 0.001) {
				     grid.text("***", x, y, gp = gpar(fontsize = 3))
			     } else if (pval < 0.01) {
				     grid.text("**", x, y, gp = gpar(fontsize = 3))
			     } else if (pval < 0.05) {
				     grid.text("*", x, y, gp = gpar(fontsize = 3))
			     }
		     },
                     name = "Z-scaled scores",
                     column_title = heatmap_title,
                     column_names_gp = gpar(fontsize = 4), 
                     row_names_gp = gpar(fontsize = 3),
                     column_names_rot = 45,
                     cluster_rows = row_dend,
                     cluster_columns = col_dend)
  
  # Draw the heatmap
  draw(heatmap, heatmap_legend_side= "bottom")
  
  return(list(heatmap_data = combined_df_heatmap, significance_tags = tags_df_heatmap))
}

heatmap_plots <- list()
for(i in 1:length(all_csv_data_1)){
  pdf(file = paste0(names(all_csv_data_1)[[i]], ".pdf"))
  heatmap_plots[[i]] <- process_and_plot_heatmap(all_csv_data_1[[i]], names(all_csv_data_1[i]))
  dev.off()
}
