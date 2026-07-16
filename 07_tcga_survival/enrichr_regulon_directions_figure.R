## Figure: top GO-BP-2023 terms per regulon direction (4 panels).
suppressMessages({library(ggplot2); library(dplyr)})
IN  <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/07_tcga_survival/output_regulon_survival/enrichr_regulon_directions"
NTOP <- 8

lab <- c(MaleUnique_Down="Male regulon (Down after KD)",
         MaleUnique_Up  ="Male regulon (Up after KD)",
         FemaleUnique_Down="Female regulon (Down after KD)",
         FemaleUnique_Up  ="Female regulon (Up after KD)")
cols <- c("Male regulon (Down after KD)"="#1E6FBF","Male regulon (Up after KD)"="#8FC1EC",
          "Female regulon (Down after KD)"="#C2185B","Female regulon (Up after KD)"="#F48FB1")

dat <- bind_rows(lapply(names(lab), function(s){
  d <- read.csv(file.path(IN, paste0("enrichr_GO_BP2023_", s, ".csv")))
  d <- d[order(d$Adjusted.P.value), ][seq_len(min(NTOP, nrow(d))), ]
  d$Term  <- sub("\\s*\\(GO:\\d+\\)$", "", d$Term)      # drop GO id
  d$Term  <- ifelse(nchar(d$Term) > 46, paste0(substr(d$Term,1,44),"…"), d$Term)
  d$Direction <- factor(lab[[s]], levels = unname(lab))
  d$neglog10 <- -log10(d$Adjusted.P.value)
  d$row <- seq_len(nrow(d))
  d
}))
## unique per-row id for y-axis (avoids duplicate term names); row1=smallest adjP -> top
dat$uid <- paste0(as.integer(dat$Direction), "_", sprintf("%02d", dat$row))
lvl <- dat %>% arrange(Direction, desc(row)) %>% pull(uid) %>% unique()
dat$uid <- factor(dat$uid, levels = lvl)
labmap <- setNames(dat$Term, as.character(dat$uid))

p <- ggplot(dat, aes(neglog10, uid, fill = Direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.6, color = "grey40") +
  facet_wrap(~Direction, scales = "free", ncol = 2) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_discrete(labels = labmap) +
  labs(x = expression(-log[10]~"(adjusted "*italic(P)*"-value)"), y = NULL,
       title = "Egr1 regulon GO enrichment by direction (GO BP 2023)") +
  theme_bw(base_size = 15) +
  theme(axis.text.y  = element_text(size = 14, color = "black"),
        axis.text.x  = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        plot.title   = element_text(size = 18, hjust = 0.5, face = "bold"),
        strip.text   = element_text(size = 15, face = "bold"),
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank())

outf <- file.path(IN, "enrichr_GO_BP2023_regulon_directions.pdf")
ggsave(outf, p, width = 15, height = 9.5, device = "pdf")
cat("wrote:", outf, "\n")
