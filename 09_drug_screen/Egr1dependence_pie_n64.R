suppressMessages(library(ggplot2))
out <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output/"
dep <- 37; ind <- 27; tot <- dep + ind          # 64
d <- data.frame(
  cat = factor(c("Egr1-dependent\n(sig. interaction)","Egr1-independent"),
               levels=c("Egr1-dependent\n(sig. interaction)","Egr1-independent")),
  n   = c(dep, ind))
d$pct <- round(100*d$n/tot)
d$lab <- paste0(d$n, "\n(", d$pct, "%)")
# position labels at slice mid-angle
d$cum <- cumsum(d$n); d$mid <- d$cum - d$n/2
p <- ggplot(d, aes(x="", y=n, fill=cat)) +
  geom_col(width=1, colour="white", linewidth=1) +
  coord_polar(theta="y", start=0, direction=1) +
  geom_text(aes(y=tot-mid, label=lab), colour="white", fontface="bold", size=7, lineheight=0.9) +
  scale_fill_manual(values=c("Egr1-dependent\n(sig. interaction)"="#1f77b4","Egr1-independent"="#b3b3b3"), name=NULL) +
  labs(title=sprintf("Egr1 dependence (n=%d)", tot)) +
  theme_void(base_size=16) +
  theme(plot.title=element_text(size=20, face="bold", hjust=0.5),
        legend.text=element_text(size=15), legend.position="right",
        legend.key.size=unit(0.9,"cm"), plot.margin=margin(10,10,10,10))
ggsave(paste0(out,"Egr1dependence_pie_n64.pdf"), p, width=7.2, height=5, device=cairo_pdf)
ggsave(paste0(out,"Egr1dependence_pie_n64.png"), p, width=7.2, height=5, dpi=200, bg="white")
cat("dependent", dep, "(", d$pct[1], "%) | independent", ind, "(", d$pct[2], "%) | n", tot, "\n")
cat("WROTE Egr1dependence_pie_n64.{pdf,png}\n")
