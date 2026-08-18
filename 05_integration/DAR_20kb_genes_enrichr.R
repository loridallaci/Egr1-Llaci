# =============================================================================
# Genes near sex-biased DARs (multiome, male vs female), + Enrichr.
#
# DARs : data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv
#        filtered at DAR_PADJ / DAR_LFC (house default 0.05 and 0.5).
# "near" = gene BODY within WINDOW bp of a DAR (maxgap on findOverlaps).
#        The house default 20,000 matches DE_near_DAR_20kb_enrichment_summary.csv.
#
# All thresholds are set in the parameter block below and every filename, plot
# title and console message is DERIVED from them, so nothing can silently keep
# saying "20kb" after the window is changed.
#
# Two flavours of gene list are produced:
#   ALL  = every gene near a male-up / female-up DAR. DARs are dense, so these
#          are close to genome-wide and Enrichr on them is uninformative;
#          off by default (RUN_ALL_SETS).
#   DEconc = sex-DE genes within WINDOW of a DIRECTION-CONCORDANT DAR
#          (male-higher genes near male-up DARs; female-higher near female-up).
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

## ---- analysis parameters (everything downstream is derived from these) --------
WINDOW   <- 20000   # bp: max distance from gene body to DAR
DAR_PADJ <- 0.05    # DAR significance
DAR_LFC  <- 0.5     # DAR |avg_log2FC|
DE_PADJ  <- 0.05    # multiome male-vs-female DE significance

## label used in every filename and plot title - derived, so it can never
## disagree with WINDOW (e.g. 20000 -> "20kb", 1500 -> "1.5kb", 500 -> "500bp")
WIN_TAG <- if (WINDOW >= 1000) {
  v <- WINDOW / 1000
  paste0(if (v == round(v)) format(v) else format(round(v, 1)), "kb")
} else paste0(WINDOW, "bp")
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
  paste0("05_integration/output/enrichment_DAR_", WIN_TAG) else
  paste0("05_integration/output/enrichment_DAR_", WIN_TAG, "_allChr"))
dir.create(OutputDirectory, showWarnings = FALSE, recursive = TRUE)
cat("AUTOSOMAL_ONLY =", AUTOSOMAL_ONLY, "-> ", OutputDirectory, "\n")

## ---- DARs -------------------------------------------------------------------
da <- read.csv(file.path(base, "data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv"),
               stringsAsFactors = FALSE)
colnames(da)[1] <- "peak"
da <- da[da$p_val_adj <= DAR_PADJ & abs(da$avg_log2FC) >= DAR_LFC, ]
cat(sprintf("DARs passing padj<=%.3g & |LFC|>=%.3g : %d", DAR_PADJ, DAR_LFC, nrow(da)),
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
cat("within", WIN_TAG, "of MALE-up DAR  :", length(genes_M_all), "genes\n")
cat("within", WIN_TAG, "of FEMALE-up DAR:", length(genes_F_all), "genes\n")

## ---- direction-concordant, DE-restricted (interpretable set) ------------------
de <- read.csv(file.path(base, "data output/DE_male_vs_female_allcells_allgenes.csv"),
               stringsAsFactors = FALSE)
de <- de[!is.na(de$gene) & !duplicated(de$gene), ]
de_sig <- de[!is.na(de$p_val_adj) & de$p_val_adj <= DE_PADJ, ]
genes_M_de <- intersect(de_sig$gene[de_sig$avg_log2FC > 0], genes_M_all)
genes_F_de <- intersect(de_sig$gene[de_sig$avg_log2FC < 0], genes_F_all)
cat("male-higher DE genes near male-up DAR    :", length(genes_M_de), "\n")
cat("female-higher DE genes near female-up DAR:", length(genes_F_de), "\n")

## ---- how many sex-DE genes sit near a DAR, by direction -----------------------
## Concordant = same direction (male-higher gene near a male-up DAR).
## Discordant = opposite direction. Reported so the concordance is auditable.
de_M <- de_sig$gene[de_sig$avg_log2FC > 0]      # male-higher
de_F <- de_sig$gene[de_sig$avg_log2FC < 0]      # female-higher
genes_any <- union(genes_M_all, genes_F_all)
de_tab <- data.frame(
  set          = c("male-higher DE", "female-higher DE"),
  n_DE         = c(length(de_M), length(de_F)),
  n_near_anyDAR    = c(sum(de_M %in% genes_any),   sum(de_F %in% genes_any)),
  n_near_concordant= c(length(genes_M_de),         length(genes_F_de)),
  n_near_discordant= c(sum(de_M %in% genes_F_all), sum(de_F %in% genes_M_all)))
de_tab$pct_near_anyDAR     <- round(100 * de_tab$n_near_anyDAR / de_tab$n_DE, 1)
de_tab$pct_near_concordant <- round(100 * de_tab$n_near_concordant / de_tab$n_DE, 1)
cat("\n-- sex-DE genes near a DAR (window", WIN_TAG, ") --\n")
print(de_tab, row.names = FALSE)
write.csv(de_tab, file.path(OutputDirectory,
          paste0("DE_near_DAR_", WIN_TAG, "_by_sex.csv")), row.names = FALSE)

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
          file.path(OutputDirectory, paste0("DAR_", WIN_TAG, "_chromosome_composition.csv")), row.names = FALSE)

## ---- write gene tables --------------------------------------------------------
wr <- function(v, f) write.csv(data.frame(SYMBOL = v), file.path(OutputDirectory, f), row.names = FALSE)
if (RUN_ALL_SETS) {
  wr(genes_M_all, paste0("genes_within", WIN_TAG, "_MALEup_DAR_all.csv"))
  wr(genes_F_all, paste0("genes_within", WIN_TAG, "_FEMALEup_DAR_all.csv"))
}
wr(genes_M_de,  paste0("genes_within", WIN_TAG, "_MALEup_DAR_DEconcordant.csv"))
wr(genes_F_de,  paste0("genes_within", WIN_TAG, "_FEMALEup_DAR_DEconcordant.csv"))

write.csv(data.frame(
  set = c("MALEup_all","FEMALEup_all","MALEup_DEconcordant","FEMALEup_DEconcordant"),
  n_genes = c(length(genes_M_all), length(genes_F_all), length(genes_M_de), length(genes_F_de)),
  n_DARs  = c(length(dar_M), length(dar_F), length(dar_M), length(dar_F)),
  window_bp = WINDOW, dar_padj = DAR_PADJ, dar_abs_log2FC = DAR_LFC,
  de_padj = DE_PADJ, autosomal_only = AUTOSOMAL_ONLY),
  file.path(OutputDirectory, paste0("DAR_", WIN_TAG, "_gene_counts.csv")), row.names = FALSE)

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
  run_enrichr(genes_M_all, paste0("MALEup_DAR_", WIN_TAG, "_all"),
              paste0("Genes within ", WIN_TAG, " of male-up DARs (all)"))
  run_enrichr(genes_F_all, paste0("FEMALEup_DAR_", WIN_TAG, "_all"),
              paste0("Genes within ", WIN_TAG, " of female-up DARs (all)"))
}
run_enrichr(genes_M_de,  paste0("MALEup_DAR_", WIN_TAG, "_DEconc"),
            paste0("Male-higher DE genes within ", WIN_TAG, " of male-up DARs"))
run_enrichr(genes_F_de,  paste0("FEMALEup_DAR_", WIN_TAG, "_DEconc"),
            paste0("Female-higher DE genes within ", WIN_TAG, " of female-up DARs"))

message("\nDone -> ", OutputDirectory)
