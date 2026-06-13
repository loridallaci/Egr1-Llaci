library(eulerr)
library(gridExtra)

# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
output_dir <- "output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_tcga_genes        <- 12701
n_sig_male_motifs   <- 567
n_sig_female_motifs <- 350
n_top30_M           <- 23
n_top30_F           <- 23
n_sig_MM            <- 23
n_sig_MF            <- 23
n_sig_FF            <- 23
n_sig_FM            <- 23

# =============================================================================
# Helper function
# =============================================================================

make_venn <- function(n_sig_motifs, n_tcga_genes, n_overlap, n_sig_multi,
                      title, motif_color = "#378ADD",
                      visual_scale = 0.05,
                      tcga_scale   = 0.01) {
  
  fit <- euler(c(
    "Sig\nMotifs"             = round(n_sig_motifs * visual_scale),
    "TCGA\nGenes"             = round(n_tcga_genes * tcga_scale),
    "Sig\nMotifs&TCGA\nGenes" = n_overlap
  ))
  
  plot(fit,
       fills      = list(fill  = c(motif_color, "#B0C4DE"), alpha = 0.7),
       edges      = list(col   = c("#1a5fa8",   "#6a8caf"), lwd   = 2),
       labels     = list(fontface = "bold", cex = 0.85),
       quantities = list(
         labels = c(
           format(n_sig_motifs - n_overlap, big.mark = ","),
           format(n_tcga_genes - n_overlap, big.mark = ","),
           paste0(n_overlap, " tested\n(", n_sig_multi, " significant)")
         ),
         cex = 0.85, fontface = "bold"
       ),
       main = list(label = title, fontface = "bold", cex = 1)
  )
}

# =============================================================================
# All 4 combined in one PDF
# =============================================================================

p_MM <- make_venn(n_sig_male_motifs,   n_tcga_genes, n_top30_M, n_sig_MM,
                  "Male motifs × Male TCGA patients",     "#378ADD",
                  visual_scale = 0.05, tcga_scale = 0.01)

p_MF <- make_venn(n_sig_male_motifs,   n_tcga_genes, n_top30_M, n_sig_MF,
                  "Male motifs × Female TCGA patients",   "#378ADD",
                  visual_scale = 0.05, tcga_scale = 0.01)

p_FF <- make_venn(n_sig_female_motifs, n_tcga_genes, n_top30_F, n_sig_FF,
                  "Female motifs × Female TCGA patients", "#D4537E",
                  visual_scale = 0.05, tcga_scale = 0.01)

p_FM <- make_venn(n_sig_female_motifs, n_tcga_genes, n_top30_F, n_sig_FM,
                  "Female motifs × Male TCGA patients",   "#D4537E",
                  visual_scale = 0.05, tcga_scale = 0.01)

cairo_pdf(file.path(output_dir, "venn_motif_tcga_4combinations.pdf"),
          width = 14, height = 10)
grid.arrange(p_MM, p_MF, p_FF, p_FM, ncol = 2)
dev.off()

# =============================================================================
# Individual PDFs
# =============================================================================

cairo_pdf(file.path(output_dir, "venn_male_motifs_male_tcga.pdf"),
          width = 7, height = 6)
print(p_MM)
dev.off()

cairo_pdf(file.path(output_dir, "venn_female_motifs_female_tcga.pdf"),
          width = 7, height = 6)
print(p_FF)
dev.off()

cat("Saved to:", output_dir, "\n")