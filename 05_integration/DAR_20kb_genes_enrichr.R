# =============================================================================
# Genes within 20 kb of sex-biased DARs (multiome, male vs female), + Enrichr.
#
# DARs      : data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv
#             house threshold p_val_adj <= 0.05 & |avg_log2FC| >= 0.5
# "within 20 kb" : gene BODY within 20,000 bp of a DAR (maxgap on findOverlaps),
#             matching the window used in DE_near_DAR_20kb_enrichment_summary.csv
# Autosomal only (chrX/chrY dropped), matching DAR_promoter_genes_by_sex.R.
#
# Two flavours of gene list are produced:
#   ALL  = every gene within 20 kb of a male-up / female-up DAR
#          (NOTE: DARs are dense, so these lists are close to genome-wide and
#           Enrichr on them is expected to be uninformative - reported anyway)
#   DEconc = sex-DE genes (padj<=0.05) within 20 kb of a DIRECTION-CONCORDANT DAR
#          (male-higher genes near male-up DARs; female-higher near female-up)
#          This is the interpretable set.
#
# Same Enrichr workflow/databases/fonts as Enrichr_female_c4_vs_big.R.
# Needs the Enrichr website (internet) - runs on the author's machine.
# =============================================================================
suppressMessages({
  library(GenomicRanges); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db); library(enrichR); library(ggplot2); library(dplyr)
})

base <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"

WINDOW <- 20000
## AUTOSOMAL_ONLY = FALSE (default) -> keep every chromosome including chrX/chrY/chrM.
##   This is the version that reproduces the canonical 1,312 / 1,219 concordant
##   gene counts, so it is the house default. Set TRUE to drop chrX/chrY/chrM.
AUTOSOMAL_ONLY <- as.logical(Sys.getenv("AUTOSOMAL_ONLY", "FALSE"))

## RUN_ALL_SETS = FALSE (default) -> only analyse the DE-concordant gene sets.
##   The "all genes within 20 kb of a DAR" sets are ~9.5k genes (about half the
##   genome, because DARs are dense) and return nothing interpretable - the
##   female set had no significant term. Not reported; off by default.
RUN_ALL_SETS <- as.logical(Sys.getenv("RUN_ALL_SETS", "FALSE"))

## REPLOT_ONLY = TRUE -> skip the Enrichr query and rebuild the figures from the
##   already-saved .txt result tables. Use when only the plot styling changed,
##   so the numbers on disk stay exactly as originally queried.
REPLOT_ONLY <- as.logical(Sys.getenv("REPLOT_ONLY", "FALSE"))

OutputDirectory <- file.path(base, if (AUTOSOMAL_ONLY)
  "05_integration/output/enrichment_DAR_20kb" else
  "05_integration/output/enrichment_DAR_20kb_allChr")
dir.create(OutputDirectory, showWarnings = FALSE, recursive = TRUE)
cat("AUTOSOMAL_ONLY =", AUTOSOMAL_ONLY, "-> ", OutputDirectory, "\n")

## ---- DARs -------------------------------------------------------------------
da <- read.csv(file.path(base, "data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv"),
               stringsAsFactors = FALSE)
colnames(da)[1] <- "peak"
da <- da[da$p_val_adj <= 0.05 & abs(da$avg_log2FC) >= 0.5, ]
cat("DARs passing padj<=0.05 & |LFC|>=0.5 :", nrow(da),
    " (male-up", sum(da$avg_log2FC > 0), "/ female-up", sum(da$avg_log2FC < 0), ")\n")

peak2gr <- function(ids) {
  m <- regmatches(ids, regexec("^(.*)-([0-9]+)-([0-9]+)$", ids))
  GRanges(vapply(m, `[`, "", 2),
          IRanges(as.integer(vapply(m, `[`, "", 3)), as.integer(vapply(m, `[`, "", 4))))
}
dar_M <- peak2gr(da$peak[da$avg_log2FC > 0])
dar_F <- peak2gr(da$peak[da$avg_log2FC < 0])

## ---- mm10 gene bodies -> symbols, autosomal only -----------------------------
gn <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene, single.strand.genes.only = TRUE)
gn <- keepStandardChromosomes(gn, pruning.mode = "coarse")   # drops unplaced scaffolds only
if (AUTOSOMAL_ONLY) gn <- gn[!as.character(seqnames(gn)) %in% c("chrX", "chrY", "chrM")]
sym <- suppressMessages(
  AnnotationDbi::select(org.Mm.eg.db, keys = gn$gene_id, keytype = "ENTREZID", columns = "SYMBOL"))
gn$SYMBOL <- sym$SYMBOL[match(gn$gene_id, sym$ENTREZID)]
gn <- gn[!is.na(gn$SYMBOL)]
cat("autosomal mm10 genes with a symbol:", length(gn), "\n")

near <- function(dar) sort(unique(gn$SYMBOL[queryHits(
  findOverlaps(gn, dar, maxgap = WINDOW, ignore.strand = TRUE))]))

genes_M_all <- near(dar_M)
genes_F_all <- near(dar_F)
cat("within 20kb of MALE-up DAR  :", length(genes_M_all), "genes\n")
cat("within 20kb of FEMALE-up DAR:", length(genes_F_all), "genes\n")

## ---- direction-concordant, DE-restricted (interpretable set) ------------------
de <- read.csv(file.path(base, "data output/DE_male_vs_female_allcells_allgenes.csv"),
               stringsAsFactors = FALSE)
de <- de[!is.na(de$gene) & !duplicated(de$gene), ]
de_sig <- de[!is.na(de$p_val_adj) & de$p_val_adj <= 0.05, ]
genes_M_de <- intersect(de_sig$gene[de_sig$avg_log2FC > 0], genes_M_all)
genes_F_de <- intersect(de_sig$gene[de_sig$avg_log2FC < 0], genes_F_all)
cat("male-higher DE genes near male-up DAR    :", length(genes_M_de), "\n")
cat("female-higher DE genes near female-up DAR:", length(genes_F_de), "\n")

## ---- chromosome composition of each gene set ----------------------------------
chrom_of <- function(v) as.character(seqnames(gn))[match(v, gn$SYMBOL)]
chrtab <- do.call(rbind, lapply(
  list(MALEup_all = genes_M_all, FEMALEup_all = genes_F_all,
       MALEup_DEconcordant = genes_M_de, FEMALEup_DEconcordant = genes_F_de),
  function(v) { ch <- chrom_of(v)
    data.frame(n_total = length(v), n_chrX = sum(ch == "chrX", na.rm = TRUE),
               n_chrY = sum(ch == "chrY", na.rm = TRUE), n_chrM = sum(ch == "chrM", na.rm = TRUE)) }))
chrtab$set <- rownames(chrtab); chrtab$pct_chrX <- round(100 * chrtab$n_chrX / chrtab$n_total, 2)
print(chrtab[, c("set","n_total","n_chrX","pct_chrX","n_chrY","n_chrM")], row.names = FALSE)
write.csv(chrtab[, c("set","n_total","n_chrX","pct_chrX","n_chrY","n_chrM")],
          file.path(OutputDirectory, "DAR_20kb_chromosome_composition.csv"), row.names = FALSE)

## ---- write gene tables --------------------------------------------------------
wr <- function(v, f) write.csv(data.frame(SYMBOL = v), file.path(OutputDirectory, f), row.names = FALSE)
if (RUN_ALL_SETS) {
  wr(genes_M_all, "genes_within20kb_MALEup_DAR_all.csv")
  wr(genes_F_all, "genes_within20kb_FEMALEup_DAR_all.csv")
}
wr(genes_M_de,  "genes_within20kb_MALEup_DAR_DEconcordant.csv")
wr(genes_F_de,  "genes_within20kb_FEMALEup_DAR_DEconcordant.csv")

write.csv(data.frame(
  set = c("MALEup_all","FEMALEup_all","MALEup_DEconcordant","FEMALEup_DEconcordant"),
  n_genes = c(length(genes_M_all), length(genes_F_all), length(genes_M_de), length(genes_F_de)),
  n_DARs  = c(length(dar_M), length(dar_F), length(dar_M), length(dar_F)),
  window_bp = WINDOW, autosomal_only = AUTOSOMAL_ONLY),
  file.path(OutputDirectory, "DAR_20kb_gene_counts.csv"), row.names = FALSE)

## ---- Enrichr ------------------------------------------------------------------
websiteLive <- getOption("enrichR.live")
if (isTRUE(websiteLive)) setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2023")   # only the db used for the manuscript figures

## ---- house-style GO bar plot -------------------------------------------------
## Horizontal bars, x = gene count, fill = adjusted p on a log10 blue->red ramp
## (blue = most significant; this direction is deliberate house style).
## Full term text is kept - the "(GO:0000000)" accession is stripped and long
## names are wrapped, so nothing is truncated.
plot_enrich_house <- function(result, title_lab, outstem) {
  d <- result[order(result$Adjusted.P.value), ]
  d <- head(d, 10)
  d$GeneCount <- lengths(strsplit(d$Genes, ";"))
  d$Term <- sub("\\s*\\(GO:[0-9]+\\)\\s*$", "", d$Term)     # drop the GO accession
  d$Term <- sub("\\s*\\(WP[0-9]+\\)\\s*$", "", d$Term)
  d$Term <- stringr::str_wrap(d$Term, width = 40)           # wrap, never truncate
  d$Term <- factor(d$Term, levels = rev(d$Term))            # most significant on top

  p <- ggplot(d, aes(x = GeneCount, y = Term, fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "red", name = "Adj. P-Value", trans = "log10") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of Genes", y = NULL, title = title_lab) +
    theme_minimal(base_size = 14) +
    theme(plot.title   = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title   = element_text(size = 16),
          axis.text.y  = element_text(size = 14, lineheight = 0.95, hjust = 1),
          axis.text.x  = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 14),
          plot.margin  = margin(12, 18, 12, 12))

  ## height scales with the number of wrapped lines so nothing looks cramped
  n_lines <- sum(vapply(strsplit(levels(d$Term), "\n"), length, integer(1)))
  h <- max(7, 0.55 * n_lines + 1.8)
  ggsave(paste0(outstem, ".pdf"), p, width = 12, height = h, bg = "white",
         device = grDevices::cairo_pdf, limitsize = FALSE)
  ggsave(paste0(outstem, ".png"), p, width = 12, height = h, dpi = 300, bg = "white",
         limitsize = FALSE)
}

run_enrichr <- function(genes, label, title_lab) {
  cat("\n== Enrichr:", label, "(", length(genes), "genes )\n")
  if (length(genes) < 5) { cat("  too few genes, skipped\n"); return(invisible(NULL)) }
  if (!REPLOT_ONLY) {
    enriched <- try(enrichr(genes, dbs), silent = TRUE)
    if (inherits(enriched, "try-error")) { cat("  ENRICHR FAILED\n"); return(invisible(NULL)) }
  }
  for (db in dbs) {
    f <- file.path(OutputDirectory, paste0(label, "_", db, ".txt"))
    if (REPLOT_ONLY) {
      if (!file.exists(f)) { cat("  no saved table for", db, "- skipped\n"); next }
      result <- read.delim(f, stringsAsFactors = FALSE)
      cat("  replotting from", basename(f), "\n")
    } else {
      result <- enriched[[db]]
    }
    if (is.null(result) || nrow(result) == 0) { cat("  skip", db, "(no results)\n"); next }
    ## row.names = FALSE: with row.names = TRUE the header has 9 fields but each
    ## data row has 10, so every column is shifted by one when the file is opened
    ## in Excel (the "Adjusted.P.value" header then sits over the raw P.value).
    if (!REPLOT_ONLY) write.table(result, f, quote = FALSE, row.names = FALSE, sep = "\t")
    top <- head(result[order(result$P.value), c("Term","Overlap","P.value","Adjusted.P.value")], 5)
    print(top, row.names = FALSE)
    plot_enrich_house(result, title_lab, file.path(OutputDirectory, paste0(label, "_", db)))
  }
}

if (RUN_ALL_SETS) {
  run_enrichr(genes_M_all, "MALEup_DAR_20kb_all",   "Genes within 20 kb of male-up DARs (all)")
  run_enrichr(genes_F_all, "FEMALEup_DAR_20kb_all", "Genes within 20 kb of female-up DARs (all)")
}
run_enrichr(genes_M_de,  "MALEup_DAR_20kb_DEconc",   "Male-higher DE genes within 20 kb of male-up DARs")
run_enrichr(genes_F_de,  "FEMALEup_DAR_20kb_DEconc", "Female-higher DE genes within 20 kb of female-up DARs")

message("\nDone -> ", OutputDirectory)
