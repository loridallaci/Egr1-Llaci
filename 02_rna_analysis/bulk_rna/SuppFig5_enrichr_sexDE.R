## =====================================================================
## enrichR GO on the Supp Fig 5c/5d sex-comparison DEGs, in EACH direction.
## Four gene sets:
##   control : higher-in-female (2,996) / higher-in-male (2,324)
##   Egr1-KD : higher-in-female (3,180) / higher-in-male (2,432)
## House style (see enrichr-house-style memory / enrichr_overlaps.R):
##   dbs = GO MF/CC/BP 2023 + WikiPathways mouse; plot GO_BP; horizontal bars,
##   x = gene count, blue(low-p)->red log10 Adjusted.P.value gradient; top 10;
##   white bg; presentation fonts. Requires internet.
## =====================================================================

suppressMessages({library(enrichR); library(ggplot2); library(dplyr); library(stringr)})
setEnrichrSite("Enrichr")

de_dir  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_sexComparison"
out_dir <- file.path(de_dir, "enrichr")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dbs     <- c("GO_Biological_Process_2023")   # only the db used for the manuscript figures
plot_db <- "GO_Biological_Process_2023"

## DE table -> direction-split gene lists (raw p<=0.05 & |lfc|>=0.5); +lfc = higher in Male
split_genes <- function(f) {
  d <- read.delim(file.path(de_dir, f), sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
  sig <- !is.na(d$pvalue) & !is.na(d$log2FoldChange) & d$pvalue <= 0.05 & abs(d$log2FoldChange) >= 0.5
  list(male   = toupper(unique(na.omit(d$SYMBOL[sig & d$log2FoldChange >=  0.5]))),
       female = toupper(unique(na.omit(d$SYMBOL[sig & d$log2FoldChange <= -0.5]))))
}
ctrl <- split_genes("Male_vs_Female_control_NoTreatg1_DE_vst_filtered.txt")
kd   <- split_genes("Male_vs_Female_Egr1KD_gRNA3_DE_vst_filtered.txt")
sets <- list(
  Control_FemaleHigher = ctrl$female, Control_MaleHigher = ctrl$male,
  Egr1KD_FemaleHigher  = kd$female,   Egr1KD_MaleHigher  = kd$male)
for (nm in names(sets)) cat(sprintf("%-22s %d genes\n", nm, length(sets[[nm]])))

run_one <- function(genes, tag) {
  genes <- genes[genes != "" & !is.na(genes)]
  res <- tryCatch(enrichr(genes, dbs), error=function(e){message("enrichr failed: ",e$message);NULL})
  if (is.null(res)) return(invisible())
  for (db in names(res))
    write.table(res[[db]], file.path(out_dir, sprintf("%s_enrichrResults_%s.txt", tag, db)),
                sep="\t", quote=FALSE, row.names=FALSE)
  plt <- res[[plot_db]] %>%
    mutate(GeneCount = str_count(Genes, ";") + 1) %>%
    arrange(Adjusted.P.value) %>% slice_head(n = 10) %>%
    mutate(Term = str_wrap(Term, width = 40), Term = factor(Term, levels = rev(Term)))
  if (nrow(plt) == 0) { cat(sprintf("  %s: no terms\n", tag)); return(invisible()) }
  ph <- min(max(7, 0.9 * nrow(plt)), 16)
  p <- ggplot(plt, aes(x = GeneCount, y = Term, fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "red", name = "Adj. P-Value", trans = "log10") +
    theme_minimal(base_size = 14) +
    labs(title = sprintf("GO Enrichment - %s", sub("_", " ", tag)), x = "Number of Genes", y = "GO Term") +
    theme(axis.text.y=element_text(size=14, hjust=1), axis.text.x=element_text(size=14),
          axis.title=element_text(size=16), legend.text=element_text(size=14),
          legend.title=element_text(size=14), plot.title=element_text(size=18, face="bold", hjust=0.5),
          plot.margin = margin(10,10,30,10))
  ggsave(file.path(out_dir, sprintf("%s_enrichrResults_%s.pdf", tag, plot_db)), p, width=10, height=ph, bg="white")
  ggsave(file.path(out_dir, sprintf("%s_enrichrResults_%s.png", tag, plot_db)), p, width=10, height=ph, dpi=200, bg="white")
  cat(sprintf("\n== %s : top GO-BP ==\n", tag))
  print(head(data.frame(Term=as.character(plt$Term), padj=signif(plt$Adjusted.P.value,3), genes=plt$Overlap), 6), row.names=FALSE)
}
for (nm in names(sets)) run_one(sets[[nm]], nm)
cat(sprintf("\nDone. enrichR results in %s\n", out_dir))
