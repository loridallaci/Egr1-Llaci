## =====================================================================
## Forest plots of RENIN-motif TF gene-expression HR in TCGA GBM, with
## THREE lines per TF (Male / Female / All patients), split into three
## panels by TF category: Shared, Male-unique, Female-unique.
## Also writes a supplementary table of the 23-per-sex TFs (categorised,
## with HR/CI/p in male, female and all TCGA patients).
## Reads cox_wide_<set>_COMBINED.csv made by motif_multivariate_TCGA.R.
## =====================================================================
suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(openxlsx)})
base   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/05_tcga_survival/motif_multivariate_TCGA"
outdir <- file.path(base, "forest_by_TFcategory"); dir.create(outdir, showWarnings = FALSE)

read_set <- function(set){
  d <- read.csv(file.path(base, set, paste0("cox_wide_", set, "_COMBINED.csv")), check.names = FALSE)
  sl <- substr(set,1,1)
  g <- function(coh, stat) suppressWarnings(as.numeric(d[[paste0("Expression_", stat, "_", sl, coh)]]))
  data.frame(TF=d$SYMBOL_UPPER,
    HR_Male=g("M","coef"),   Lo_Male=g("M","Lower95"),   Hi_Male=g("M","Upper95"),   p_Male=g("M","pvalue"),
    HR_Female=g("F","coef"), Lo_Female=g("F","Lower95"), Hi_Female=g("F","Upper95"), p_Female=g("F","pvalue"),
    HR_All=g("A","coef"),    Lo_All=g("A","Lower95"),    Hi_All=g("A","Upper95"),    p_All=g("A","pvalue"),
    ovp_Male=as.numeric(d$Overall_pvalue_Males), ovp_Female=as.numeric(d$Overall_pvalue_Females),
    stringsAsFactors=FALSE)
}
M <- read_set("MaleMotifs"); Fm <- read_set("FemaleMotifs")
shared <- intersect(M$TF, Fm$TF); male_only <- setdiff(M$TF, Fm$TF); female_only <- setdiff(Fm$TF, M$TF)

dat <- bind_rows(M %>% filter(TF %in% c(shared, male_only)), Fm %>% filter(TF %in% female_only))
dat <- dat[!duplicated(dat$TF), ]
dat$Category <- ifelse(dat$TF %in% shared, "Shared",
               ifelse(dat$TF %in% male_only, "Male-unique", "Female-unique"))
dat$ordp <- ifelse(dat$Category=="Female-unique", dat$ovp_Female, dat$ovp_Male)

long <- bind_rows(
  transmute(dat, TF, Category, ordp, Cohort="Male",   HR=HR_Male,   Lo=Lo_Male,   Hi=Hi_Male,   p=p_Male),
  transmute(dat, TF, Category, ordp, Cohort="Female", HR=HR_Female, Lo=Lo_Female, Hi=Hi_Female, p=p_Female),
  transmute(dat, TF, Category, ordp, Cohort="All",    HR=HR_All,    Lo=Lo_All,    Hi=Hi_All,    p=p_All))
long$Cohort <- factor(long$Cohort, levels=c("All","Female","Male"))  # Male dodged on TOP

make_plot <- function(cat, fname, title){
  d <- long %>% filter(Category==cat, is.finite(HR), is.finite(Lo), is.finite(Hi))
  o <- d %>% distinct(TF, ordp) %>% arrange(desc(ordp))               # smallest p on top
  d$TF  <- factor(d$TF, levels=o$TF)
  d$sig <- ifelse(!is.na(d$p) & d$p<0.05, "*", "")
  p <- ggplot(d, aes(HR, TF, color=Cohort)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=Lo, xmax=Hi), height=0.3, position=position_dodge(width=0.7)) +
    geom_point(size=2.8, position=position_dodge(width=0.7)) +
    geom_text(aes(x=Hi, label=sig), position=position_dodge(width=0.7), hjust=-0.3, vjust=0.75, size=6, show.legend=FALSE) +
    scale_color_manual(values=c("Male"="#1E90FF","Female"="#FF69B4","All"="#33A02C"), name=NULL) +
    scale_x_log10() + scale_y_discrete(expand=expansion(add=0.7)) +
    guides(color=guide_legend(reverse=TRUE)) +                        # legend: Male, Female, All
    labs(x="RENIN motif TF gene-expression HR (95% CI, log scale)", y="RENIN motif TF",
         title=title, subtitle="Male / Female / All TCGA GBM patients (multivariate Cox); * = HR p < 0.05") +
    theme_bw() + theme(legend.position="top")
  ggsave(file.path(outdir, fname), p, width=8, height=max(3, 0.42*length(unique(d$TF))+1.6))
  cat(sprintf("%-14s %2d TFs -> %s\n", cat, length(unique(d$TF)), fname))
}
make_plot("Shared",        "forest_SharedTFs_MFA.pdf",    "Shared RENIN motif TFs - TCGA GBM")
make_plot("Male-unique",   "forest_MaleUnique_MFA.pdf",   "Male-unique RENIN motif TFs - TCGA GBM")
make_plot("Female-unique", "forest_FemaleUnique_MFA.pdf", "Female-unique RENIN motif TFs - TCGA GBM")

## ---- supplementary table (all 39 TFs; 23 per sex = 7 shared + 16 unique) ----
supp <- dat %>% transmute(Category, TF,
  HR_male=round(HR_Male,3),   CI_low_male=round(Lo_Male,3),   CI_high_male=round(Hi_Male,3),   p_male=signif(p_Male,3),
  HR_female=round(HR_Female,3),CI_low_female=round(Lo_Female,3),CI_high_female=round(Hi_Female,3),p_female=signif(p_Female,3),
  HR_all=round(HR_All,3),     CI_low_all=round(Lo_All,3),     CI_high_all=round(Hi_All,3),     p_all=signif(p_All,3)) %>%
  arrange(factor(Category, levels=c("Shared","Male-unique","Female-unique")), p_male)
write.csv(supp, file.path(outdir,"SupplTable_TCGA_TFs_shared_unique.csv"), row.names=FALSE)
openxlsx::write.xlsx(supp, file.path(outdir,"SupplTable_TCGA_TFs_shared_unique.xlsx"))
cat("\nShared:", length(shared), " Male-unique:", length(male_only), " Female-unique:", length(female_only),
    " | total", nrow(dat), "\nWrote plots + SupplTable to:", outdir, "\n")
