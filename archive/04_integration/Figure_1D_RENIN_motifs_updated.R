# =============================================================================
# 03_atac_analysis: Motif Enrichment Analysis - Female and Male CREs
# =============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(ggrepel)

# --- Load object and CRE lists ------------------------------------------------

obj         <- readRDS("/home/lllaci/data/obj_RENIN_processed.rds")
cre_lists   <- readRDS("/home/lllaci/data/cre_lists_female_male_RENIN.rds")
female_cres <- cre_lists$female_cres
male_cres   <- cre_lists$male_cres

DefaultAssay(obj) <- "peaks"

# =============================================================================
# Helper: motif enrichment plot
# =============================================================================

plot_motif_enrichment <- function(motifs, title, num_top_label = 30,
                                  sig_color = "dodgerblue3",
                                  nonsig_color = "grey80",
                                  xlim_fixed = NULL,
                                  ylim_fixed = NULL,
                                  cutoff_line = NULL,
                                  cutoff_color = "black") {
  
  motifs$logp        <- -log10(motifs$pvalue)
  motifs$significant <- motifs$p.adjust <= 0.05
  
  motifs$color <- ifelse(
    motifs$significant,
    sig_color,
    nonsig_color
  )
  
  motifs$border <- ifelse(motifs$significant, 0.4, 0.1)
  
  motifs$label <- ifelse(
    rank(-motifs$logp) <= num_top_label,
    motifs$motif.name,
    ""
  )
  
  max_y <- if (!is.null(ylim_fixed)) ylim_fixed else max(motifs$logp, na.rm = TRUE) * 1.05
  max_x <- if (!is.null(xlim_fixed)) xlim_fixed else max(motifs$fold.enrichment, na.rm = TRUE) * 1.05
  
  p <- ggplot(motifs, aes(x = fold.enrichment, y = logp, label = label)) +
    geom_point(
      aes(fill = color),
      color  = "black",
      stroke = motifs$border,
      pch    = 21,
      size   = 2
    ) +
    scale_fill_identity() +
    
    # --- ADD CUT-OFF LINE ---
    {
      if (!is.null(cutoff_line)) {
        geom_hline(yintercept = cutoff_line,
                   linetype = "dashed",
                   color = cutoff_color)
      }
    } +
    
    geom_text_repel(
      max.overlaps       = 500,
      size               = 3,
      point.padding      = 0.5,
      force              = 2,
      box.padding        = 0.4,
      min.segment.length = 0,
      segment.color      = "grey40",
      segment.size       = 0.4,
      segment.alpha      = 0.8
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = title,
      x     = "Fold enrichment",
      y     = expression(-log[10](p-value))
    ) +
    ylim(c(0, max_y)) +
    xlim(c(0, max_x)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text  = element_text(color = "black", size = 14)
    )
  
  return(p)
}

# =============================================================================
# MALE-ENRICHED CREs
# =============================================================================

message("Running motif enrichment for male-enriched CREs...")
set.seed(12345)

motifs_male <- FindMotifs(object = obj, features = male_cres)
motifs_male$logp <- -log10(motifs_male$pvalue)
motifs_male <- motifs_male[order(motifs_male$logp, decreasing = TRUE), ]

# --- cutoff for top 30 motifs (male) ---
cutoff_male <- motifs_male$logp[30]

cat("EGR1 rank in male motifs:", which(motifs_male$motif.name == "EGR1"), "\n")
print(motifs_male[motifs_male$motif.name == "EGR1", ])

cat("Significant male motifs (p.adjust <= 0.05):",
    sum(motifs_male$p.adjust <= 0.05), "\n")

write.csv(motifs_male, "M_all_motifs_updated.csv", row.names = TRUE)

cairo_pdf("motif_enrichment_mal_updatede.pdf", width = 7, height = 6)
print(
  plot_motif_enrichment(
    motifs_male,
    title = "Motif Enrichment: Male-enriched CREs",
    sig_color = "dodgerblue3",
    nonsig_color = "grey85",
    cutoff_line = cutoff_male,
    cutoff_color = "red"
  )
)
dev.off()

message("Done. Male motif enrichment plot saved.")

# =============================================================================
# FEMALE-ENRICHED CREs
# =============================================================================

message("Running motif enrichment for female-enriched CREs...")
set.seed(12345)

motifs_female <- FindMotifs(object = obj, features = female_cres)
motifs_female$logp <- -log10(motifs_female$pvalue)
motifs_female <- motifs_female[order(motifs_female$logp, decreasing = TRUE), ]

# --- cutoff for top 30 motifs (female) ---
cutoff_female <- motifs_female$logp[30]

cat("EGR1 rank in female motifs:", which(motifs_female$motif.name == "EGR1"), "\n")
print(motifs_female[motifs_female$motif.name == "EGR1", ])

cat("Significant female motifs (p.adjust <= 0.05):",
    sum(motifs_female$p.adjust <= 0.05), "\n")

write.csv(motifs_female, "F_all_motifs_updated.csv", row.names = TRUE)

cairo_pdf("motif_enrichment_female_updated.pdf", width = 7, height = 6)
print(
  plot_motif_enrichment(
    motifs_female,
    title = "Motif Enrichment: Female-enriched CREs",
    sig_color = "pink",
    nonsig_color = "grey85",
    cutoff_line = cutoff_female,
    cutoff_color = "blue"
  )
)
dev.off()

message("Done. Female motif enrichment plot saved.")