# run_all.R  — source from data_analysis_visualization/ or use Rscript run_all.R

if (interactive()) {
  tryCatch(
    setwd(dirname(rstudioapi::getSourceEditorContext()$path)),
    error = function(e) message("Set working directory to data_analysis_visualization/ manually.")
  )
}

library(tidyverse)
library(ggpubr)

source("load_data.R")

old_figs <- list.files(FIG_DIR, full.names = TRUE)
if (length(old_figs) > 0) file.remove(old_figs)

source("analysis.R")
