# Data directory structure

library(here)

data_directory <- here("data")
output_directory <- here("output")
figure_directory <-  here("figures")
de_directory <- here("diffexp")
r_script_directory <- here("R")

for(directory in c(data_directory, output_directory, figure_directory,de_directory, r_script_directory)){
  if (!dir.exists(directory)) { dir.create(directory) }
}

# Default values for Differential Expression Testing
max_FDR <- 0.05
min_LFC <- 1

# Libraries

# Graphics
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(plotly)
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
palette(pal_nejm()(8))
palette_cb_gray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette_cb_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Formated tables
library(DT)
library(htmltools)

# annotation
library(biomaRt)

# Differential Expression
library(limma)
library(edgeR)
library(Glimma)

# heatmaps
library(pheatmap)

# PCA plots
library(genefilter)

# helper functions
source(file.path(r_script_directory,"01-build_gene_annotation.R"))
source(file.path(r_script_directory,"02-qc_sequence_effort.R"))
source(file.path(r_script_directory,"03-filter_lowly_expressed.R"))
source(file.path(r_script_directory,"04-PCA_functions.R"))
source(file.path(r_script_directory,"05-volcano_plot.R"))
source(file.path(r_script_directory,"06-Glimma_functions.R"))
source(file.path(r_script_directory,"07-TF_Ligand_Receptor_functions.R"))
