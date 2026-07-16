## =====================================================================
## Supplementary Fig 1b - RNA UMAP from the RPCA-INTEGRATED multiome object.
##
## Plots the STORED umap.RPCA embedding carried in the processed lot6 object
## (the same embedding used for panel 1e, the EGR1 motif panel). Reading the
## stored embedding - rather than recomputing UMAP - is deterministic, so
## panels 1b and 1e sit on an identical layout cell-for-cell. No coordinate
## files: the embedding comes straight out of the Seurat object.
##
## Single panel, male + female overlaid, coloured by sex.
## Runs locally (Seurat attaches fine for read-only; the Matrix/chromVAR
## segfault only affects RunChromVAR, which we do not call here).
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(ggplot2) })

rds <- "C:/Users/loril/Downloads/lot6_merged_RPCAintegrated_011524_motifsAdded_030124_chromVARadded_030124_V4_050424_2x.rds"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/multiome_rna/output/figures/suppl_fig1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sex_cols <- c(female = "#F39AC9", male = "#4A6FE3")

obj <- readRDS(rds)
## sex coded 1 = female, 2 = male (Cell Ranger ARC -1/-2 aggregation order)
obj$sex <- factor(ifelse(obj$sex %in% c("1", 1), "female", "male"), levels = c("female","male"))

emb <- Embeddings(obj, "umap.RPCA")
df  <- data.frame(umapRPCA_1 = emb[,1], umapRPCA_2 = emb[,2], sex = obj$sex)
## randomize draw order so neither sex is plotted entirely on top of the other
set.seed(42); df <- df[sample(nrow(df)), ]

p <- ggplot(df, aes(umapRPCA_1, umapRPCA_2, colour = sex)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_colour_manual(values = sex_cols, name = "Sex") +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(title = "RNA UMAP (RPCA-integrated)", x = "umapRPCA_1", y = "umapRPCA_2") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14))
ggsave(file.path(out, "SupplFig1b_RNA_UMAP_RPCA.pdf"), p, width = 6.5, height = 5.5)

cat("Done. RNA UMAP (single, by sex) ->", out, "\n")
