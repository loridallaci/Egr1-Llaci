# MAKE DAR VOLCANO ANALYSIS
suppressMessages(library(Seurat)) 
suppressMessages(library(Signac))
suppressMessages(library(SeuratWrappers)) 
suppressMessages(library(RENIN)) 
suppressMessages(library(harmony))
#plan(sequential)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)

# original (author's machine): "/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds"
lot6 = readRDS("output/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds")
obj <- lot6
DefaultAssay(obj) <- 'RNA'

# Normalize 
obj <- NormalizeData(obj)
# Find variable features
obj <- FindVariableFeatures(obj)
# Scale data (important!)
obj <- ScaleData(obj)

# PCA -> neighbors -> clusters -> UMAP
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

# Add sex info
obj$sex <- ifelse(grepl("-1$", colnames(obj)), "female", 
                  ifelse(grepl("-2$", colnames(obj)), "male", NA))

# Process the PEAK data
DefaultAssay(obj) <- "peaks"
obj<- FindTopFeatures(obj, min.cutoff = 5)
obj<- RunTFIDF(obj)
obj<- RunSVD(obj)


obj<- RunUMAP(object = obj, reduction = 'lsi', dims = 2:20)
obj<- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:20)
obj<- FindClusters(object = obj, verbose = FALSE, algorithm = 3)

Idents(obj) <- 'sex'
DefaultAssay(obj) <- 'peaks'

# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = obj,
  ident.1 = "male",
  ident.2 = "female",
  test.use = 'wilcox',
  min.pct = 0.1, 
  logfc.threshold = 0
)

head(da_peaks)
#write.csv(da_peaks, "lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv", row.names = TRUE)
write.csv(da_peaks, "lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks_logfc0.csv", row.names = TRUE)




























#Made the volcano in R - In Multiome_update_01292025.Rmd file

11/10/25

# DO everything again, but with DAR on normalized ATACseq
# I already run DAR from peaks, M vs F. Now find genes nearest those peaks, and then overlap with DE genes. 
# Let's plot DAR volcano first
```{r}
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
da_peaks <- read.csv('output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv')  # original: 'C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Code for figures/Files associated with code/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks.csv'


make_gr <- function(df) {
  parts <- str_split_fixed(da_peaks$X, "-", 3)
  GRanges(
    seqnames = parts[,1],
    ranges = IRanges(
      start = as.numeric(parts[,2]),
      end = as.numeric(parts[,3])
    )
  )
}

da_peaks_gr <- make_gr(da_peaks)
da_peaks_gr

# Remove "chr" prefix from da_gr
seqlevels(da_peaks_gr) <- sub("^chr", "", seqlevels(da_peaks_gr))
# Check again
seqlevels(da_peaks_gr)[1:5]
# Should now match genes_gr: "1" "10" "11" "12" "13"


# Get genes from the annotation
edb <- EnsDb.Mmusculus.v79
genes_gr <- genes(edb)

# Find nearest genes
nearest_genes <- nearest(da_peaks_gr, genes_gr)
nearest_genes_total <- genes_gr[nearest_genes]
nearest_genes_total
# To know the distane to nearest gene
nearest_genes_distance <- distanceToNearest(da_peaks_gr, genes_gr)
nearest_genes_distance
```
# add gene name to original peaks df
```{r}
nearest_gene_names <- mcols(nearest_genes_total)$gene_name
da_peaks$nearest_gene <- nearest_gene_names
da_peaks$nearest_gene_id <- mcols(nearest_genes_total)$gene_id
da_peaks
```

# The original code had   da_peaks$avg_log2FC >= 0.05 & da_peaks$p_val_adj <= 0.05, I'll change the av log2FC to > 0.5

# make volcano plot
```{r}
library(dplyr)
da_peaks$significance <- "Not Sig"
da_peaks$significance[
  da_peaks$avg_log2FC >= 0.5 & da_peaks$p_val_adj <= 0.05
] <- "Male-up"
da_peaks$significance[
  da_peaks$avg_log2FC <= -0.5 & da_peaks$p_val_adj <= 0.05
] <- "Female-up"

library(ggrepel)

p <- ggplot(da_peaks, aes(
  x = avg_log2FC,
  y = -log10(p_val_adj),
  color = significance
)) +
  geom_point(
    alpha = 0.7,
    size = 1.8,
    shape = 16
  ) +
  scale_color_manual(values = c(
    "Male-up" = "blue",
    "Female-up" = "pink",
    "Not Sig" = "grey70"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "Differentially Accessible Peaks",
    x = "log2 Fold Change",
    y = "-log10(adj. p-value)"
  )
p

ggsave(
  filename = "output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks_volcano_log2fc05_padj005_010726.png",  # original: "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Figure1_Jan2026/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks_volcano_log2fc05_padj005_010726.png"
  plot = p,
  width = 7,
  height = 6,
  dpi = 600,
  units = "in"
)

```

# Count them
```{r}
n_male_sig <- sum(da_peaks$significance == "Male-up")
n_female_sig <- sum(da_peaks$significance == "Female-up")

n_total_sig <- n_male_sig + n_female_sig

n_male_sig
n_female_sig
n_total_sig
table(da_peaks$significance)


```
# save da_peaks to implort into RENIN terminal
```{r}
# Save
write.csv(da_peaks, "output/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks_edited.csv", row.names = FALSE)  # original: "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Code for figures/Files associated with code/lot6_MvsF_PeaksNormalized_DifferentialyAccesible_peaks_edited.csv"

```