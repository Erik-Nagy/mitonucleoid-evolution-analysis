# run_all.R  — source from data_analysis_visualization/ or use Rscript run_all.R

library(tidyverse)
library(ggpubr)

source("00_load_data.R")

old_figs <- list.files(FIG_DIR, full.names = TRUE)
if (length(old_figs) > 0) file.remove(old_figs)

source("analysis.R")
