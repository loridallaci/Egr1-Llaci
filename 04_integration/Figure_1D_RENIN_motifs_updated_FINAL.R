# =============================================================================
# 03_atac_analysis: Motif Enrichment Analysis - Female and Male CREs
# =============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(ggrepel)

# --- Load object and CRE lists ------------------------------------------------

obj         <- readRDS("output/obj_RENIN_processed.rds")  # original: "/home/lllaci/data/obj_RENIN_processed.rds"
cre_lists   <- readRDS("output/cre_lists_female_male_RENIN.rds")  # original: "/home/lllaci/data/cre_lists_female_male_RENIN.rds"
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
      segment.alpha      = 0.8,
      ylim               = c(cutoff_line, NA) # labels cannot go below cutoff line
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

cairo_pdf("motif_enrichment_male_updated.pdf", width = 7, height = 6)
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
    cutoff_color = "red"
  )
)
dev.off()

message("Done. Female motif enrichment plot saved.")


## EXTRA
# =============================================================================
# Overlap of significant motifs between male and female CREs
# =============================================================================
# 316 motifs were shared between the two datasets. 
sig_male   <- motifs_male  [motifs_male$p.adjust   <= 0.05, ]
sig_female <- motifs_female[motifs_female$p.adjust <= 0.05, ]

shared_ids <- intersect(sig_male$motif, sig_female$motif)
male_only  <- setdiff(sig_male$motif,   sig_female$motif)
female_only <- setdiff(sig_female$motif, sig_male$motif)

cat("Significant motifs in male CREs:    ", nrow(sig_male),    "\n")
cat("Significant motifs in female CREs:  ", nrow(sig_female),  "\n")
cat("Shared (significant in both):       ", length(shared_ids),  "\n")
cat("Male-only (significant only in M):  ", length(male_only),   "\n")
cat("Female-only (significant only in F):", length(female_only), "\n")

# Show the shared motif names
shared_table <- merge(
  sig_male[,   c("motif", "motif.name", "fold.enrichment", "p.adjust")],
  sig_female[, c("motif", "motif.name", "fold.enrichment", "p.adjust")],
  by = "motif", suffixes = c("_male", "_female")
)
shared_table <- shared_table[order(shared_table$p.adjust_male), ]

print(head(shared_table, 20))

write.csv(shared_table, "shared_significant_motifs_M_vs_F.csv", row.names = FALSE)
# --- Unique motifs per sex ----------------------------------------------------

male_only_table <- sig_male[sig_male$motif %in% male_only,
                            c("motif", "motif.name", "fold.enrichment",
                              "pvalue", "p.adjust")]
male_only_table <- male_only_table[order(male_only_table$p.adjust), ]

female_only_table <- sig_female[sig_female$motif %in% female_only,
                                c("motif", "motif.name", "fold.enrichment",
                                  "pvalue", "p.adjust")]
female_only_table <- female_only_table[order(female_only_table$p.adjust), ]

cat("\nTop 10 male-only significant motifs:\n")
print(head(male_only_table, 10))

cat("\nTop 10 female-only significant motifs:\n")
print(head(female_only_table, 10))

write.csv(male_only_table,
          "male_only_significant_motifs.csv",   row.names = FALSE)
write.csv(female_only_table,
          "female_only_significant_motifs.csv", row.names = FALSE)

## top 30
# --- Overlap of top 30 motifs (by p-value rank) ------------------------------

top_n <- 30

# motifs_male and motifs_female are already sorted by logp descending
top30_male   <- head(motifs_male,   top_n)
top30_female <- head(motifs_female, top_n)

shared_top30   <- intersect(top30_male$motif, top30_female$motif)
top30_male_only   <- setdiff(top30_male$motif,   top30_female$motif)
top30_female_only <- setdiff(top30_female$motif, top30_male$motif)

cat("\nTop", top_n, "motif overlap:\n")
cat("  Shared in both top", top_n, ":     ", length(shared_top30),     "\n")
cat("  Top", top_n, "male-only:           ", length(top30_male_only),   "\n")
cat("  Top", top_n, "female-only:         ", length(top30_female_only), "\n")

# Build a side-by-side table of the shared top-30 motifs
shared_top30_table <- merge(
  top30_male[,   c("motif", "motif.name", "fold.enrichment",
                   "pvalue", "p.adjust", "logp")],
  top30_female[, c("motif", "motif.name", "fold.enrichment",
                   "pvalue", "p.adjust", "logp")],
  by = "motif", suffixes = c("_male", "_female")
)
shared_top30_table <- shared_top30_table[order(shared_top30_table$logp_male,
                                               decreasing = TRUE), ]

cat("\nShared top-30 motifs:\n")
print(shared_top30_table[, c("motif", "motif.name_male",
                             "fold.enrichment_male", "p.adjust_male",
                             "fold.enrichment_female", "p.adjust_female")])

write.csv(shared_top30_table,
          "shared_top30_motifs_M_vs_F.csv", row.names = FALSE)

# Also save the sex-specific top-30 lists for completeness
write.csv(top30_male[top30_male$motif %in% top30_male_only, ],
          "top30_male_only_motifs.csv", row.names = FALSE)
write.csv(top30_female[top30_female$motif %in% top30_female_only, ],
          "top30_female_only_motifs.csv", row.names = FALSE)

