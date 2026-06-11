library(Seurat)

# Load object
lot6 <- readRDS("/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds")
obj <- lot6
# Make sure RNA is normalized
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)

# Set identity to sex
Idents(obj) <- "sex"

# Run DE: female vs male
de_sex <- FindMarkers(
  object    = obj,
  ident.1   = "male",
  ident.2   = "female",
  assay     = "RNA",
  test.use  = "wilcox",    # standard for scRNA-seq
  min.pct   = 0.1,         # gene must be detected in at least 10% of cells
  logfc.threshold = 0   # minimum log fold change
)

# Add gene names as column
de_sex$gene <- rownames(de_sex)

# Save ALL genes (unfiltered)
write.csv(de_sex, "DE_male_vs_female_allcells_allgenes.csv", row.names = FALSE)

# Filter for significant genes
de_sig <- de_sex[de_sex$p_val_adj <= 0.05, ]

# Split into female-higher and male-higher
male_higher   <- de_sig[de_sig$avg_log2FC >=  0.5, ]
female_higher <- de_sig[de_sig$avg_log2FC <= -0.5, ]

# Save results
write.csv(de_sig,        "DE_male_vs_female_allcells_significant.csv",  row.names = FALSE)
write.csv(male_higher,   "DE_male_vs_female_allcells_maleHigher.csv",   row.names = FALSE)
write.csv(female_higher, "DE_male_vs_female_allcells_femaleHigher.csv", row.names = FALSE)

# Quick summary
cat("Total genes tested:      ", nrow(de_sex),     "\n")
cat("Total significant:       ", nrow(de_sig),      "\n")
cat("Higher in male:          ", nrow(male_higher),   "\n")
cat("Higher in female:        ", nrow(female_higher), "\n")
