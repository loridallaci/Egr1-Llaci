# Supp Fig 2a — cortical nuclei by stage x sex, 100% paper cells (39,069).
# Counts read directly from the cortex_development_paperCells objects (.181).
suppressMessages(library(ggplot2))
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/06_cortex_development/output"

df <- data.frame(
  stage = rep(c("LateFetal","Infant","Child","Adol","Adult"), each=2),
  sex   = rep(c("male","female"), 5),
  Freq  = c(3102,4754, 3683,2723, 5418,3693, 3188,5112, 5855,1541))
df$stage <- factor(df$stage, levels=c("LateFetal","Infant","Child","Adol","Adult"))
df$sex   <- factor(df$sex,   levels=c("male","female"))
stopifnot(sum(df$Freq)==39069)

sex_cols <- c(male="#4A6FE3", female="#F39AC9")

pc <- ggplot(df, aes(stage, Freq, fill=sex)) +
  geom_col(position=position_dodge(0.75), width=0.7, colour="black") +
  geom_text(aes(label=Freq), position=position_dodge(0.75), vjust=-0.4, size=4.6) +
  scale_fill_manual(values=sex_cols, name="Sex", labels=c("Male","Female")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.10))) +
  labs(title="Cortical nuclei by developmental stage and sex", x=NULL, y="Number of nuclei") +
  theme_classic(base_size=16) +
  theme(plot.title   = element_text(size=18, face="bold"),
        axis.title   = element_text(size=16),
        axis.text    = element_text(size=14),
        axis.text.x  = element_text(angle=15, hjust=1),
        legend.title = element_text(size=14),
        legend.text  = element_text(size=14),
        legend.position="top")

ggsave(file.path(out,"SupplFig2_cortex_cellnumbers.pdf"), pc, width=8, height=5.5)
ggsave(file.path(out,"SupplFig2_cortex_cellnumbers.png"), pc, width=8, height=5.5, dpi=150)
cat("wrote SupplFig2_cortex_cellnumbers.pdf/png; total =", sum(df$Freq), "\n")
