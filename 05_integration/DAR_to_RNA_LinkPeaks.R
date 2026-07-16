## =====================================================================
## Sex-biased DAR peaks -> genes by CO-ACCESSIBILITY (Signac LinkPeaks) ->
## concordance with RNA sex-bias.
##
## Companion to DAR_to_RNA_nearestGene.R. Nearest-gene linkage is unfair to
## DISTAL (female-biased) peaks, which regulate genes they are not nearest to.
## LinkPeaks correlates each peak's accessibility with nearby gene expression
## (<=500 kb) across cells, so distal peaks get linked to their actual targets.
##
## Convention: positive avg_log2FC = MALE-higher in both DAR and RNA DE tables.
##
## Inputs:
##   the processed multiome object (peaks assay has Annotation + GC stats)
##   data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv
##   data output/DE_male_vs_female_allcells_allgenes.csv
## Outputs (05_integration/output/):
##   DAR_to_RNA_LinkPeaks_links.csv        (all significant DAR peak-gene links)
##   DAR_to_RNA_LinkPeaks_summary.csv      (concordance, nearest-gene vs LinkPeaks)
##   Fig_DAR_to_RNA_LinkPeaks_concordance.pdf
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(GenomicRanges); library(ggplot2)
                   library(ChIPseeker); library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db) })
base <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
rds  <- "C:/Users/loril/Downloads/lot6_merged_RPCAintegrated_011524_motifsAdded_030124_chromVARadded_030124_V4_050424_2x.rds"
out  <- file.path(base, "05_integration/output"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

da <- read.csv(file.path(base, "data output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv"), stringsAsFactors = FALSE)
de <- read.csv(file.path(base, "data output/DE_male_vs_female_allcells_allgenes.csv"), stringsAsFactors = FALSE)
de <- de[!is.na(de$gene) & !duplicated(de$gene), ]; rownames(de) <- de$gene

## ---- classify DAR peaks: male-up PROMOTER vs female-up DISTAL -------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
g <- function(ids){ m <- regmatches(ids, regexec("^(.*)-([0-9]+)-([0-9]+)$", ids))
  GRanges(vapply(m,`[`,"",2), IRanges(as.integer(vapply(m,`[`,"",3)), as.integer(vapply(m,`[`,"",4))), peakid = ids) }
mu <- da[[1]][da$avg_log2FC >=  0.5 & da$p_val_adj <= 0.05]
fu <- da[[1]][da$avg_log2FC <= -0.5 & da$p_val_adj <= 0.05]
am <- as.data.frame(annotatePeak(g(mu), TxDb = txdb, tssRegion = c(-2000,500), verbose = FALSE))
af <- as.data.frame(annotatePeak(g(fu), TxDb = txdb, tssRegion = c(-2000,500), verbose = FALSE))
male_prom_peaks   <- am$peakid[grepl("^Promoter", am$annotation)]
female_dist_peaks <- af$peakid[grepl("^Distal Intergenic", af$annotation)]
cat(sprintf("male-up promoter peaks=%d ; female-up distal peaks=%d\n", length(male_prom_peaks), length(female_dist_peaks)))

## ---- object + candidate genes within 500 kb of those DARs ----------
## NB: the DAR table is a DIFFERENT (normalized) peak set than the object's
## `peaks` assay, so DAR peaks are matched to object peaks by OVERLAP, not name.
obj <- readRDS(rds); DefaultAssay(obj) <- "peaks"
male_prom_gr   <- granges(g(male_prom_peaks))
female_dist_gr <- granges(g(female_dist_peaks))
darGR <- c(male_prom_gr, female_dist_gr)

ann <- Annotation(obj[["peaks"]])
genes_gr <- range(split(ann, ann$gene_name))                 # gene-body span per gene
genes_gr <- unlist(genes_gr); genes_gr$gene_name <- names(genes_gr)
near <- unique(subsetByOverlaps(genes_gr, resize(darGR, width(darGR) + 1e6, fix = "center"))$gene_name)
genes.use <- intersect(near, rownames(obj[["RNA"]]))
cat("candidate genes within 500kb of DARs, expressed:", length(genes.use), "\n")

## ---- LinkPeaks -----------------------------------------------------
DefaultAssay(obj) <- "peaks"
obj <- LinkPeaks(obj, peak.assay = "peaks", expression.assay = "RNA",
                 genes.use = genes.use, distance = 5e5, min.cells = 10)
lk <- as.data.frame(Links(obj[["peaks"]]))
lk <- lk[lk$pvalue <= 0.05, ]
cat("significant links (p<=0.05):", nrow(lk), "\n")

## keep DAR-peak links; tag DAR class by OVERLAP (DAR peakset != object peakset)
psep <- if (any(grepl(":", rownames(obj[["peaks"]])))) c(":","-") else c("-","-")
lkgr <- StringToGRanges(lk$peak, sep = psep)
lk$class <- ifelse(overlapsAny(lkgr, male_prom_gr),    "Male-up PROMOTER",
             ifelse(overlapsAny(lkgr, female_dist_gr), "Female-up DISTAL", NA))
lk <- lk[!is.na(lk$class), ]
lk$rna_lfc <- de[lk$gene, "avg_log2FC"]; lk$rna_padj <- de[lk$gene, "p_val_adj"]
lk <- lk[!is.na(lk$rna_lfc), ]
write.csv(lk, file.path(out, "DAR_to_RNA_LinkPeaks_links.csv"), row.names = FALSE)

## ---- concordance per class (unique linked genes) -------------------
bg_male <- mean(de$avg_log2FC > 0)
summ <- do.call(rbind, lapply(split(lk, lk$class), function(s){
  gg <- unique(s$gene); sub <- de[gg, ]
  a <- sum(sub$avg_log2FC > 0); b <- nrow(sub) - a
  cM <- sum(de$avg_log2FC > 0) - a; d <- nrow(de) - a - b - cM
  data.frame(class = s$class[1], n_links = nrow(s), n_genes = nrow(sub),
             pct_male_hi = 100*mean(sub$avg_log2FC > 0),
             pct_female_hi = 100*mean(sub$avg_log2FC < 0),
             fisher_male_p   = fisher.test(matrix(c(a,b,cM,d),2), alternative="greater")$p.value,
             fisher_female_p = fisher.test(matrix(c(b,a,d,cM),2), alternative="greater")$p.value) }))
write.csv(summ, file.path(out, "DAR_to_RNA_LinkPeaks_summary.csv"), row.names = FALSE)
cat(sprintf("\nBackground male-higher = %.1f%%\n", 100*bg_male)); print(summ)

## ---- figure: concordant fraction, LinkPeaks --------------------------
plt <- rbind(data.frame(class = summ$class, dir = "male-higher",   pct = summ$pct_male_hi),
             data.frame(class = summ$class, dir = "female-higher", pct = summ$pct_female_hi))
p <- ggplot(plt, aes(class, pct, fill = dir)) +
  geom_col(position = position_dodge(0.75), width = 0.7, colour = "black") +
  geom_hline(yintercept = 100*bg_male, linetype = "dashed") +
  annotate("text", x = 0.6, y = 100*bg_male + 3, label = "background\nmale-higher", size = 4, hjust = 0) +
  scale_fill_manual(values = c("male-higher" = "#4A6FE3", "female-higher" = "#F39AC9"), name = NULL) +
  labs(title = "DAR co-accessibility links (LinkPeaks): RNA sex-bias",
       x = NULL, y = "% of linked genes") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.position = "top")
ggsave(file.path(out, "Fig_DAR_to_RNA_LinkPeaks_concordance.pdf"), p, width = 8, height = 5.5)
cat("Done ->", out, "\n")
