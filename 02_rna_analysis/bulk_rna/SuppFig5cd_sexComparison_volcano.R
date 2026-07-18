## =====================================================================
## Supp Fig 5c/5d : Male-vs-Female bulk-RNA differential expression, in
##   control (WT) cells and in Egr1-KD cells.
##
## Uses the SAME DESeq2 pipeline as DE_Egr1KD_gRNA3_vs_Neg1.R (design ~0+group,
## same outlier/Neg2 removal), then extracts the two SEX contrasts (rather than
## the within-sex KD contrasts) and writes their DE tables + volcano plots.
##
##   5c  control : NoTreatment_gRNA1  Male vs Female  -> 2,996 F-higher / 2,324 M-higher
##   5d  Egr1-KD : Egr1KD_gRNA3       Male vs Female  -> 3,180 F-higher / 2,432 M-higher
##   (verified 2026-07-18; raw p <= 0.05 & |log2FC| >= 0.5)
##
## Volcano styling matches Egr1KD_volcano.R (log2FC capped at +/-4 for display,
## classify on real values, label sex-marker genes + Egr1).
## =====================================================================

suppressMessages({library(DESeq2); library(dplyr); library(ggplot2); library(ggrepel)})

dd  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/data_raw"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_sexComparison"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

## ---- same setup as DE_Egr1KD_gRNA3_vs_Neg1.R -----------------------------
countData  <- read.delim(file.path(dd,"Dedup_Counts.txt"), header=TRUE, stringsAsFactors=FALSE, row.names=1)
geneData   <- read.delim(file.path(dd,"GeneInfo_Egr1KD.csv"), header=TRUE, stringsAsFactors=FALSE, row.names=1, sep=',')
geneData   <- geneData[!duplicated(geneData$ensembl_gene_id), ]
sampleData <- read.delim(file.path(dd,"Meta_data_Egr1KD.csv"), header=TRUE, row.names=1, sep=',')
sampleData <- as.data.frame(t(sampleData)); colnames(sampleData)[3] <- 'group'
countData  <- countData[, rownames(sampleData)]

samples_to_remove <- c("F6_Egr1_gRNA2_Rep2","F6_Egr1_gRNA2_Rep3","F6_Egr1_gRNA2_Rep4")
sampleData <- sampleData %>% filter(!(rownames(.) %in% samples_to_remove)) %>%
  filter(Treatment != "NoTreatment_gRNA2")
countData  <- countData[, colnames(countData) %in% rownames(sampleData)]
geneData   <- geneData[geneData$ensembl_gene_id %in% rownames(countData), ]
sampleData[,c('Treatment','Sex','group')] <- lapply(sampleData[,c('Treatment','Sex','group')], factor)

dds <- DESeqDataSetFromMatrix(countData, sampleData, design = ~0 + group)
rownames(geneData) <- geneData$ensembl_gene_id
mcols(dds) <- DataFrame(mcols(dds), geneData[rownames(dds), ])
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

## ---- volcano (styling from Egr1KD_volcano.R) -----------------------------
LOG2FC_CAP <- 4
label_genes <- c("Egr1","Xist","Ddx3y","Uty","Eif2s3y","Kdm5d")   # Egr1 + canonical sex markers

make_de_volcano <- function(a, b, tag, title, out_txt, out_pdf) {
  res <- results(dds, contrast = c("group", a, b))               # +lfc = higher in a (Male)
  res$SYMBOL <- mcols(dds)$external_gene_name
  write.table(as.data.frame(res), file.path(out, out_txt), quote=FALSE, sep="\t", col.names=NA)

  de <- as.data.frame(res); de <- de[!is.na(de$pvalue), ]
  de$pvalue[de$pvalue == 0] <- 1e-300
  de$lfc_plot <- pmax(pmin(de$log2FoldChange, LOG2FC_CAP), -LOG2FC_CAP)
  up   <- "Higher in Male";   dn <- "Higher in Female"
  de$cls <- "NO"
  de$cls[de$log2FoldChange >=  0.5 & de$pvalue <= 0.05] <- up
  de$cls[de$log2FoldChange <= -0.5 & de$pvalue <= 0.05] <- dn
  n_up <- sum(de$cls==up); n_dn <- sum(de$cls==dn)
  de$lab <- ifelse(de$SYMBOL %in% label_genes, de$SYMBOL, NA)
  cols <- setNames(c("#1E90FF","#FF69B4","grey60"), c(up,dn,"NO"))

  p <- ggplot(de, aes(lfc_plot, -log10(pvalue), col=cls, label=lab)) +
    geom_point(alpha=0.8, shape=16) +
    geom_text_repel(box.padding=2.5, max.overlaps=Inf, size=4, fontface="italic", show.legend=FALSE) +
    scale_color_manual(values=cols, name="Expression Change") +
    scale_x_continuous(limits=c(-LOG2FC_CAP,LOG2FC_CAP), breaks=seq(-LOG2FC_CAP,LOG2FC_CAP,1)) +
    theme_minimal(base_size=14) +
    labs(title=title,
         subtitle=sprintf("%s: %d  |  %s: %d  |  log2FC capped at ±%d", dn, n_dn, up, n_up, LOG2FC_CAP),
         x="Log2 Fold Change (Male / Female)", y="-Log10 P-value") +
    theme(plot.title=element_text(hjust=0.5, size=16, face="bold"),
          plot.subtitle=element_text(hjust=0.5, size=11, colour="grey40"),
          axis.title=element_text(size=14), axis.text=element_text(size=12),
          legend.text=element_text(size=12), legend.title=element_text(size=13))
  cairo_pdf(file.path(out, out_pdf), width=8, height=6); print(p); dev.off()
  ggsave(file.path(out, sub("\\.pdf$",".png",out_pdf)), p, width=8, height=6, dpi=200, bg="white")
  cat(sprintf("%-22s female-higher %d  male-higher %d\n", tag, n_dn, n_up))
}

make_de_volcano("NoTreatment_gRNA1_Male","NoTreatment_gRNA1_Female","SuppFig5c control M-vs-F",
  "Male vs Female - control (WT) cells",
  "Male_vs_Female_control_NoTreatg1_DE_vst_filtered.txt",
  "SuppFig5c_MalevsFemale_control_volcano.pdf")

make_de_volcano("Egr1KD_gRNA3_Male","Egr1KD_gRNA3_Female","SuppFig5d KD M-vs-F",
  "Male vs Female - Egr1 knockdown cells",
  "Male_vs_Female_Egr1KD_gRNA3_DE_vst_filtered.txt",
  "SuppFig5d_MalevsFemale_Egr1KD_volcano.pdf")

cat("\nDone. DE tables + volcanos in:\n  ", out, "\n")
