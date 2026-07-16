## =====================================================================
## Sex-biased DAR peaks -> nearest gene -> concordance with RNA sex-bias.
##
## Tests whether accessibility differences translate to expression differences,
## separately for male-up PROMOTER peaks and female-up DISTAL/intergenic peaks.
## Peaks are linked to genes by NEAREST-GENE annotation (ChIPseeker); a proper
## co-accessibility version is in DAR_to_RNA_LinkPeaks.R.
##
## Convention: positive avg_log2FC = MALE-higher (ident.1 = male) in BOTH the
## DAR table and the RNA DE table.
##
## Inputs (data output/):
##   lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv   (DARs)
##   DE_male_vs_female_allcells_allgenes.csv                      (RNA DE)
## Outputs (05_integration/output/):
##   Fig_DAR_to_RNA_concordance.pdf
##   DAR_to_RNA_nearestGene_summary.csv
##
## Seurat-free (ChIPseeker/GenomicRanges only) -> runs locally.
## =====================================================================
suppressMessages({ library(ChIPseeker); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
                   library(org.Mm.eg.db); library(GenomicRanges); library(ggplot2) })
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
base <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
out  <- file.path(base, "05_integration/output"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

da <- read.csv(file.path(base, "data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv"), stringsAsFactors = FALSE)
de <- read.csv(file.path(base, "data output/DE_male_vs_female_allcells_allgenes.csv"), stringsAsFactors = FALSE)
de <- de[!is.na(de$gene) & !duplicated(de$gene), ]; rownames(de) <- de$gene

## peak id "chr-start-end" -> GRanges
g <- function(ids){ m <- regmatches(ids, regexec("^(.*)-([0-9]+)-([0-9]+)$", ids))
  GRanges(vapply(m, `[`, "", 2), IRanges(as.integer(vapply(m, `[`, "", 3)), as.integer(vapply(m, `[`, "", 4)))) }
anno <- function(ids) as.data.frame(annotatePeak(g(ids), TxDb = txdb, tssRegion = c(-2000, 500),
                                                 annoDb = "org.Mm.eg.db", verbose = FALSE))
mu <- da[[1]][da$avg_log2FC >=  0.5 & da$p_val_adj <= 0.05]   # male-up DARs
fu <- da[[1]][da$avg_log2FC <= -0.5 & da$p_val_adj <= 0.05]   # female-up DARs
am <- anno(mu); af <- anno(fu)
genes_of <- function(df, pat) unique(na.omit(df$SYMBOL[grepl(pat, df$annotation)]))
sets <- list("Male-up ALL"           = unique(na.omit(am$SYMBOL)),
             "Male-up PROMOTER"      = genes_of(am, "^Promoter"),
             "Female-up ALL"         = unique(na.omit(af$SYMBOL)),
             "Female-up DISTAL"      = genes_of(af, "^Distal Intergenic"))

bg <- mean(de$avg_log2FC > 0)
rows <- lapply(names(sets), function(nm){
  gg <- sets[[nm]]; sub <- de[gg[gg %in% de$gene], ]
  a <- sum(sub$avg_log2FC > 0); b <- nrow(sub) - a
  c <- sum(de$avg_log2FC > 0) - a; d <- nrow(de) - a - b - c
  data.frame(set = nm, n = nrow(sub),
             pct_male_hi   = 100 * mean(sub$avg_log2FC > 0),
             pct_sig_male  = 100 * mean(sub$avg_log2FC > 0 & sub$p_val_adj <= 0.05),
             pct_sig_female= 100 * mean(sub$avg_log2FC < 0 & sub$p_val_adj <= 0.05),
             fisher_male_p = fisher.test(matrix(c(a,b,c,d),2), alternative = "greater")$p.value) })
summ <- do.call(rbind, rows)
write.csv(summ, file.path(out, "DAR_to_RNA_nearestGene_summary.csv"), row.names = FALSE)
print(summ)

## ---- figure: % of linked genes significantly male-/female-higher ---
plt <- rbind(data.frame(set = summ$set, dir = "sig. male-higher",   pct = summ$pct_sig_male),
             data.frame(set = summ$set, dir = "sig. female-higher", pct = summ$pct_sig_female))
plt$set <- factor(plt$set, levels = names(sets))
p <- ggplot(plt, aes(set, pct, fill = dir)) +
  geom_col(position = position_dodge(0.75), width = 0.7, colour = "black") +
  scale_fill_manual(values = c("sig. male-higher" = "#4A6FE3", "sig. female-higher" = "#F39AC9"), name = NULL) +
  labs(title = "DAR-linked genes: RNA sex-bias concordance",
       x = NULL, y = "% linked genes DE in RNA (padj ≤ 0.05)") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"), axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 20, hjust = 1), axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14), legend.position = "top")
ggsave(file.path(out, "Fig_DAR_to_RNA_concordance.pdf"), p, width = 8, height = 5.5)
cat(sprintf("\nBackground male-higher = %.1f%%\nDone -> %s\n", 100*bg, out))
