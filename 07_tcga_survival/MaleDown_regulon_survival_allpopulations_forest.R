## Male induced (Down-after-KD) regulon survival across TCGA populations.
## Uses ONLY the committed summary values (raw-p regulon; nothing recomputed).
suppressMessages({library(ggplot2); library(readr)})
R  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci"
s  <- suppressWarnings(read_csv(file.path(R,"07_tcga_survival/output_regulon_survival/TCGA_regulon_UNIQUE_UpDown_zscore_survival_summary.csv"), show_col_types=FALSE))
d  <- as.data.frame(s[s$Geneset=="MaleUnique_Down", c("Dataset","N_genes","HR","Lower95","Upper95","Pvalue")])

nm <- c(All="All patients", Males="Male patients", Females="Female patients")
d$Dataset <- factor(nm[d$Dataset], levels = c("All patients","Female patients","Male patients"))
d$HR <- as.numeric(d$HR); d$Lower95 <- as.numeric(d$Lower95); d$Upper95 <- as.numeric(d$Upper95)
fmtp <- function(p) ifelse(p<1e-3, sprintf("%.1e", p), sprintf("%.3f", p))
d$lab <- sprintf("HR = %.2f (%.2f–%.2f),  p = %s", d$HR, d$Lower95, d$Upper95, fmtp(d$Pvalue))
xmax <- max(d$Upper95)*1.05

p <- ggplot(d, aes(HR, Dataset)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.14, linewidth = 0.9, color = "#4A4A4A") +
  geom_point(size = 6, shape = 18, color = "#B2182B") +
  geom_text(aes(x = xmax, label = lab), hjust = 1, vjust = -1.1, size = 5) +
  scale_x_continuous(trans = "log10", breaks = c(0.5,1,2,4,8), limits = c(0.5, xmax)) +
  labs(x = "Hazard ratio (95% CI)", y = NULL,
       title = "Male induced regulon (Down after KD)",
       subtitle = "High expression predicts worse survival in all TCGA-GBM populations") +
  theme_bw(base_size = 15) +
  theme(axis.text.y  = element_text(size = 16, color = "black"),
        axis.text.x  = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        plot.title   = element_text(size = 18, face = "bold"),
        plot.subtitle= element_text(size = 14, color = "grey30"),
        panel.grid.minor = element_blank())

outf <- file.path(R,"07_tcga_survival/output_regulon_survival/MaleDown_regulon_survival_allpopulations_forest.pdf")
ggsave(outf, p, width = 9, height = 4.6, device = "pdf")
cat("wrote:", basename(outf), "\n"); print(d[,c("Dataset","N_genes","HR","Lower95","Upper95","Pvalue")])
