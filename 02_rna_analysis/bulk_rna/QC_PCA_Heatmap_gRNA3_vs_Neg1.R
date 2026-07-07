## PCA + Heatmap ONLY (no DE) of the bulk RNA samples used for gRNA3-vs-Neg1.
## Removes ALL gRNA2 (Egr1KD_gRNA2 + NoTreatment_gRNA2) and the manual outliers,
## leaving NoTreatment_gRNA1 + Egr1KD_gRNA3 (both sexes). vst-transformed.
suppressMessages({library(DESeq2); library(pheatmap); library(ggplot2)})

dd  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/data_raw"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_QC_PCA_Heatmap"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cnt <- read.delim(file.path(dd, "Dedup_Counts.txt"), row.names = 1, check.names = FALSE)
sd  <- as.data.frame(t(read.delim(file.path(dd, "Meta_data_Egr1KD.csv"), row.names = 1, sep = ",")))
colnames(sd)[3] <- "group"

## remove manual outliers + all gRNA2 (Egr1 KD gRNA2 and NoTreatment gRNA2)
outliers <- c("F6_Egr1_gRNA2_Rep2", "F6_Egr1_gRNA2_Rep3", "F6_Egr1_gRNA2_Rep4")
sd <- sd[!(rownames(sd) %in% outliers), ]
sd <- sd[sd$Treatment %in% c("NoTreatment_gRNA1", "Egr1KD_gRNA3"), ]   # drops all gRNA2
cnt <- cnt[, rownames(sd)]
sd[c("Treatment","Sex","group")] <- lapply(sd[c("Treatment","Sex","group")], factor)
cat("Samples left:", nrow(sd), "\n"); print(table(sd$group))

d <- DESeqDataSetFromMatrix(round(as.matrix(cnt)), sd, ~1)
d <- d[rowSums(counts(d)) > 1, ]
v <- vst(d, blind = TRUE)

## colors
group_palette <- c(Egr1KD_gRNA3_Female = "#B266FF", Egr1KD_gRNA3_Male = "#FFB84D",
                   NoTreatment_gRNA1_Male = "#1E90FF", NoTreatment_gRNA1_Female = "#FF69B4")
sex_palette <- c(Female = "#F08080", Male = "darkgreen")

## ---- Heatmap (sample-sample Euclidean distance on vst) ----
dists <- dist(t(assay(v)))
ann <- as.data.frame(colData(v)[, c("group","Sex")])
ann$group <- droplevels(factor(ann$group)); ann$Sex <- droplevels(factor(ann$Sex))
pdf(file.path(out, "Heatmap_gRNA3_vs_Neg1.pdf"), width = 11, height = 9)
pheatmap(as.matrix(dists), col = colorRampPalette(c("blue","white","red"))(20),
         border_color = "black", annotation_col = ann,
         annotation_colors = list(Sex = sex_palette, group = group_palette),
         clustering_distance_rows = dists, clustering_distance_cols = dists,
         show_rownames = TRUE, fontsize_row = 7, fontsize_col = 7,
         main = "Sample distance (vst) - Egr1KD gRNA3 + Neg gRNA1")
invisible(dev.off())

## ---- PCA (colored by group; top 2000 variable genes) ----
pd <- plotPCA(v, intgroup = "group", ntop = 2000, returnData = TRUE)
pd$group <- factor(pd$group, levels = names(group_palette))
pv <- round(100 * attr(pd, "percentVar"))
p <- ggplot(pd, aes(PC1, PC2, fill = group)) +
  geom_point(shape = 21, size = 4, color = "black") +
  geom_text(aes(label = name), size = 2.2, vjust = -1.1) +
  scale_fill_manual(values = group_palette) +
  xlab(paste0("PC1: ", pv[1], "%")) + ylab(paste0("PC2: ", pv[2], "%")) +
  ggtitle("PCA (vst) - Egr1KD gRNA3 + Neg gRNA1") + theme_bw()
pdf(file.path(out, "PCA_gRNA3_vs_Neg1.pdf"), width = 9, height = 7); print(p); invisible(dev.off())

cat("\nPCA + Heatmap written to:\n  ", out, "\n")
