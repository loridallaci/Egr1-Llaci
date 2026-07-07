## Ortholog (babelgene) matching of RENIN top-30 motif TFs to TCGA genes,
## compared to the original uppercase method. Mirrors 01_load_and_prepare_tcga_data_updated.R
## (top 30 by p.adjust -> motif.name -> strip parens -> split "::").
suppressMessages({library(dplyr); library(tidyr); library(babelgene)})

md  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/07_tcga_survival/data_motifs"
exp <- "C:/Users/loril/Documents/multivariate_analysis/glioVis/2024-06-04_TCGA_GBM_expression.txt"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/07_tcga_survival/output_RENIN_motifs_ortholog"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

tcga <- toupper(trimws(gsub('"', '', strsplit(readLines(exp, 1), "\t")[[1]][-1])))

get_tfs <- function(f) {
  read.csv(file.path(md, f)) %>%
    arrange(p.adjust) %>% slice(1:30) %>%
    mutate(TF_name = trimws(gsub("\\s*\\([^)]*\\)", "", motif.name))) %>%
    separate_rows(TF_name, sep = "::") %>% pull(TF_name) %>% trimws() %>% unique()
}

ortho_match <- function(tf) {
  o   <- tryCatch(orthologs(genes = tf, species = "mouse", human = FALSE), error = function(e) NULL)
  hum <- if (!is.null(o) && nrow(o)) toupper(o$human_symbol) else character(0)
  mp  <- if (!is.null(o) && nrow(o)) o$symbol               else character(0)
  intersect(unique(c(hum, toupper(setdiff(tf, mp)))), tcga)
}

for (s in c("Male", "Female")) {
  tf <- get_tfs(if (s == "Male") "M_all_motifs_updated.csv" else "F_all_motifs_updated.csv")
  up <- intersect(toupper(tf), tcga)
  o  <- ortho_match(tf)
  cat(sprintf("%s: uppercase=%d | ortholog=%d | gained=%d | lost=%d\n",
              s, length(up), length(o), length(setdiff(o, up)), length(setdiff(up, o))))
  if (length(setdiff(o, up))) cat("   gained:", paste(setdiff(o, up), collapse = ", "), "\n")
  if (length(setdiff(up, o))) cat("   lost:",   paste(setdiff(up, o), collapse = ", "), "\n")
  writeLines(sort(o), file.path(out, paste0(s, "_motif_TFs_in_TCGA_ortholog.txt")))
}
