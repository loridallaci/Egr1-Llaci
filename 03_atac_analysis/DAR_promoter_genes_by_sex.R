## =====================================================================
## Genes whose PROMOTERS are sex-biased DAR peaks (male-up vs female-up).
## Promoter = ChIPseeker annotation "Promoter" (TSS -2000/+500).
## Output: one gene table per sex, with the peak's accessibility log2FC and
## (where available) the gene's male-vs-female RNA log2FC for context.
## Seurat-free -> runs locally.
## =====================================================================
suppressMessages({ library(ChIPseeker); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
                   library(org.Mm.eg.db) })
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
base <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
out  <- file.path(base, "03_atac_analysis/output/tables"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

da <- read.csv(file.path(base, "data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv"), stringsAsFactors = FALSE)
colnames(da)[1] <- "peak"
de <- read.csv(file.path(base, "data output/DE_male_vs_female_allcells_allgenes.csv"), stringsAsFactors = FALSE)
de <- de[!is.na(de$gene) & !duplicated(de$gene), ]; rownames(de) <- de$gene

g <- function(ids){ m <- regmatches(ids, regexec("^(.*)-([0-9]+)-([0-9]+)$", ids))
  GenomicRanges::GRanges(vapply(m,`[`,"",2), IRanges::IRanges(as.integer(vapply(m,`[`,"",3)), as.integer(vapply(m,`[`,"",4))), peak = ids) }

promoter_genes <- function(peakids, tag){
  a <- as.data.frame(annotatePeak(g(peakids), TxDb = txdb, tssRegion = c(-2000,500),
                                  annoDb = "org.Mm.eg.db", verbose = FALSE))
  a <- a[grepl("^Promoter", a$annotation), ]
  a <- a[!(as.character(a$seqnames) %in% c("chrX","chrY")), ]   # drop sex-chromosome genes (autosomal only)
  a$peak_log2FC <- da$avg_log2FC[match(a$peak, da$peak)]
  a$rna_log2FC  <- de[a$SYMBOL, "avg_log2FC"]
  a$rna_padj    <- de[a$SYMBOL, "p_val_adj"]
  a <- a[!is.na(a$SYMBOL), c("SYMBOL","peak","annotation","distanceToTSS","peak_log2FC","rna_log2FC","rna_padj")]
  a <- a[order(-abs(a$peak_log2FC)), ]
  a <- a[!duplicated(a$SYMBOL), ]                      # one row per gene (strongest promoter peak)
  write.csv(a, file.path(out, paste0("DAR_promoter_genes_autosomal_", tag, ".csv")), row.names = FALSE)
  cat(sprintf("\n=== %s: %d autosomal promoter genes -> DAR_promoter_genes_autosomal_%s.csv ===\n", tag, nrow(a), tag))
  cat("top 30 by |peak log2FC|:\n"); cat(paste(head(a$SYMBOL,30), collapse=", "), "\n")
  a
}

mu <- da$peak[da$avg_log2FC >=  0.5 & da$p_val_adj <= 0.05]
fu <- da$peak[da$avg_log2FC <= -0.5 & da$p_val_adj <= 0.05]
am <- promoter_genes(mu, "MALE_up")
af <- promoter_genes(fu, "FEMALE_up")

## overlap / sex-specific
cat("\nshared promoter genes (both sexes):", length(intersect(am$SYMBOL, af$SYMBOL)), "\n")
cat("male-only promoter genes:  ", length(setdiff(am$SYMBOL, af$SYMBOL)), "\n")
cat("female-only promoter genes:", length(setdiff(af$SYMBOL, am$SYMBOL)), "\n")
cat("\nDone ->", out, "\n")
