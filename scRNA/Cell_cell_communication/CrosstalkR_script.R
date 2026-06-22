library(CrossTalkeR)
library(igraph)
library(stringr)

paths <- c(
  "Condition1" = "/path/to/condition1/EV_LR.csv",
  "Condition2" = "/path/to/condition2/ThPO_LR.csv"
)

output_path <- "/path/to/output/folder/"
data <- generate_report(paths,
                out_path=output,
                out_file = 'vignettes_example.html',
                output_fmt = "html_document",
                report = TRUE,
                org = "mmu", 
                filtered_net=TRUE)