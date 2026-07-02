## =====================================================================
## RENIN motif TFs -> TCGA GBM multivariate Cox survival
## Replicates the github-code MaleMotifs/ and FemaleMotifs/ output folders.
##
## Gene sets = top-30 RENIN motifs per sex -> motif.name -> strip parens ->
##   split "::" -> unique TF names -> KEEP only those present in TCGA (31 -> 23).
## For each motif set (MaleMotifs, FemaleMotifs) x cohort (Males, Females, All):
##   per-gene multivariate Cox  Surv(survival,status) ~ Expression +
##     Recurrence + Age + Subtype + MGMT_status
## Outputs per set (mirrors the original folder layout):
##   cox_long_<set>_<cohort>.csv          (Gene, Variable, HR, Lower95, Upper95, SE, Pvalue)
##   cox_wide_<set>_<cohort>.csv          (1 row/TF; suffix MM/MF/MA, FM/FF/FA)
##   cox_wide_<set>_COMBINED.csv          (3 cohorts joined)
##   cox_wide_<set>_COMBINED_allTFs.xlsx  (pretty: Motif + Expression HR/CI/SE/p, male & female)
##   forest_plots/forest_<TF>_<cohort>.pdf
##   km_plots/KM_<TF>_<cohort>.pdf
## =====================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(survival); library(ggplot2)
})
have_surv <- requireNamespace("survminer", quietly = TRUE)
have_xlsx <- requireNamespace("openxlsx",  quietly = TRUE)

## ---- Paths ----------------------------------------------------------
base   <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/05_tcga_survival"
md     <- file.path(base, "data_motifs")
outtop <- file.path(base, "motif_multivariate_TCGA")
dir.create(outtop, showWarnings = FALSE, recursive = TRUE)
glio_dir     <- "C:/Users/loril/Documents/multivariate_analysis/glioVis"
tcga_rds_dir <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/github code/output"

all_covarID <- c("Recurrence", "Age", "Subtype", "MGMT_status")
OS_var <- "survival"; OS_event <- "status"

## ---- 1. TCGA: All / male / female pheno (RDS fast-path) --------------
fact <- function(df){ for(c in c("Recurrence","Subtype","MGMT_status")) df[[c]] <- factor(df[[c]]); df }
rds_a <- file.path(tcga_rds_dir, "tcga_pheno_maleANDfemale_allgenes_updated.rds")
rds_m <- file.path(tcga_rds_dir, "tcga_pheno_male_allgenes_updated.rds")
rds_f <- file.path(tcga_rds_dir, "tcga_pheno_female_allgenes_updated.rds")
if (all(file.exists(c(rds_a, rds_m, rds_f)))) {
  tcga_all <- fact(readRDS(rds_a)); tcga_m <- fact(readRDS(rds_m)); tcga_f <- fact(readRDS(rds_f))
} else {
  pheno <- read.table(file.path(glio_dir,"2024-06-04_TCGA_GBM_pheno.txt"), sep="\t", header=TRUE, check.names=FALSE)
  rownames(pheno) <- pheno$Sample
  ex <- read.table(file.path(glio_dir,"2024-06-04_TCGA_GBM_expression.txt"), sep="\t", header=TRUE, check.names=FALSE)
  rownames(ex) <- ex$Sample; ex <- as.matrix(ex[, setdiff(colnames(ex),"Sample")])
  common <- intersect(pheno$Sample, rownames(ex)); pheno <- pheno[common,]; ex <- ex[common,]
  pheno <- pheno[pheno$IDH1_status=="Wild-type" & !is.na(pheno$Gender),]; ex <- ex[pheno$Sample,]
  tcga_all <- fact(cbind(pheno, as.data.frame(ex)))
  tcga_m <- tcga_all[tcga_all$Gender=="Male",]; tcga_f <- tcga_all[tcga_all$Gender=="Female",]
}
cat("TCGA:", nrow(tcga_all), "All |", nrow(tcga_m), "M |", nrow(tcga_f), "F\n")

tcga_cols     <- colnames(tcga_all)
upper_to_orig <- setNames(tcga_cols, toupper(tcga_cols)); upper_to_orig <- upper_to_orig[!duplicated(names(upper_to_orig))]

## ---- 2. motif TFs per sex (top30 -> split :: -> unique -> in TCGA) ---
read_motif_tfs <- function(f) {
  d <- read.csv(file.path(md, f)) %>% arrange(p.adjust) %>% slice(1:30) %>%
    mutate(TF = trimws(gsub("\\s*\\([^)]*\\)", "", motif.name))) %>%
    separate_rows(TF, sep = "::") %>% mutate(TF = trimws(TF))
  d
}
to_tcga <- function(tf_upper) unique(unname(upper_to_orig[tf_upper[tf_upper %in% names(upper_to_orig)]]))

motif_tbl <- list(MaleMotifs   = read_motif_tfs("M_all_motifs_updated.csv"),
                  FemaleMotifs = read_motif_tfs("F_all_motifs_updated.csv"))
genesets  <- lapply(motif_tbl, function(d) to_tcga(toupper(unique(d$TF))))
for (s in names(genesets)) cat(s, ": ", length(unique(motif_tbl[[s]]$TF)), " unique TFs -> ",
                               length(genesets[[s]]), " in TCGA\n", sep = "")

## TF -> source motif.name lookup (first motif containing it), for the pretty table
tf2motif <- function(set) {
  d <- motif_tbl[[set]]
  d %>% mutate(TF_up = toupper(TF)) %>% distinct(TF_up, .keep_all = TRUE) %>%
    dplyr::select(TF_up, motif.name)
}

## ---- 3. add median-split _binary cols (for KM) ----------------------
add_binary <- function(df, genes) {
  for (g in genes) {
    if (!g %in% colnames(df)) next
    b <- paste0(g, "_binary"); if (b %in% colnames(df)) next
    df[[b]] <- factor(ifelse(df[[g]] >= median(df[[g]], na.rm = TRUE), "High", "Low"), levels = c("Low","High"))
  }
  df
}
allg <- unique(unlist(genesets))
tcga_all <- add_binary(tcga_all, allg); tcga_m <- add_binary(tcga_m, allg); tcga_f <- add_binary(tcga_f, allg)

## ---- 4. Cox helpers (same logic as Regulon_TCGA_multivariate_052726.R)
relabel_long <- function(v) {
  v <- gsub("MGMT_statusUnmethylated", "MGMT (unmethylated)",    v)
  v <- gsub("RecurrenceRecurrent",     "Recurrence (recurrent)", v)
  v <- gsub("RecurrenceSecondary",     "Recurrence (secondary)", v)
  v <- gsub("SubtypeMesenchymal",      "Mesenchymal Subtype",    v)
  v <- gsub("SubtypeProneural",        "Proneural Subtype",      v); v
}
relabel_wide <- function(v, gene) {
  v <- gsub(gene, "Expression", v, fixed = TRUE)
  v <- gsub("MGMT_statusUnmethylated", "MGMT_Unmethylated",    v)
  v <- gsub("RecurrenceRecurrent",     "Recurrence_Recurrent", v)
  v <- gsub("RecurrenceSecondary",     "Recurrence_Secondary", v)
  v <- gsub("SubtypeMesenchymal",      "Subtype_Mesenchymal",  v)
  v <- gsub("SubtypeProneural",        "Subtype_Proneural",    v); v
}
cox_long <- function(pheno, genes) {
  res <- list()
  for (g in genes) {
    if (!g %in% colnames(pheno)) next
    fm <- as.formula(paste0("Surv(",OS_var,",",OS_event,") ~ `", g, "` + ", paste(all_covarID, collapse="+")))
    fit <- tryCatch(coxph(fm, data = pheno), error = function(e) NULL); if (is.null(fit)) next
    s <- summary(fit); co <- as.data.frame(s$coefficients); ci <- as.data.frame(s$conf.int)
    for (v in rownames(co)) {
      lab <- if (v == paste0("`",g,"`") || v == g) g else relabel_long(v)
      res[[paste(g,v)]] <- data.frame(Gene=g, Variable=lab,
        HR=ci[v,"exp(coef)"], Lower95=ci[v,"lower .95"], Upper95=ci[v,"upper .95"],
        SE=co[v,"se(coef)"], Pvalue=co[v,"Pr(>|z|)"], stringsAsFactors=FALSE)
    }
  }
  bind_rows(res)
}
cox_wide <- function(pheno, genes, suffix) {
  wide <- list()
  for (g in genes) {
    if (!g %in% colnames(pheno)) next
    fm <- as.formula(paste0("Surv(",OS_var,",",OS_event,") ~ `", g, "` + ", paste(all_covarID, collapse="+")))
    fit <- tryCatch(coxph(fm, data = pheno), error = function(e) NULL); if (is.null(fit)) next
    s <- summary(fit); co <- as.data.frame(s$coefficients); ci <- as.data.frame(s$conf.int)
    row <- data.frame(SYMBOL_UPPER=g, Overall_pvalue=as.numeric(s$logtest["pvalue"]), stringsAsFactors=FALSE)
    for (v in rownames(co)) {
      lb <- relabel_wide(gsub("`","",v), g)
      row[[paste0(lb,"_coef_",suffix)]]    <- ci[v,"exp(coef)"]
      row[[paste0(lb,"_Lower95_",suffix)]] <- ci[v,"lower .95"]
      row[[paste0(lb,"_Upper95_",suffix)]] <- ci[v,"upper .95"]
      row[[paste0(lb,"_SE_",suffix)]]      <- co[v,"se(coef)"]
      row[[paste0(lb,"_pvalue_",suffix)]]  <- co[v,"Pr(>|z|)"]
    }
    wide[[g]] <- row
  }
  bind_rows(wide)
}
forest_plot <- function(gene, long, ttl) {
  gd <- long %>% filter(Gene == gene); if (!nrow(gd)) return(NULL)
  ord <- c(gene,"Age","MGMT (unmethylated)","Recurrence (recurrent)","Recurrence (secondary)",
           "Mesenchymal Subtype","Proneural Subtype")
  ord <- ord[ord %in% gd$Variable]
  gd <- gd %>% mutate(sig=ifelse(Pvalue<0.05,"p<0.05","ns"), Variable=factor(Variable, levels=rev(ord)))
  ggplot(gd, aes(HR, Variable)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40") +
    geom_point(aes(color=sig), size=4) +
    geom_errorbarh(aes(xmin=Lower95, xmax=Upper95), height=0.3) +
    scale_color_manual(values=c("p<0.05"="red","ns"="black")) +
    theme_classic(base_size=14) +
    theme(legend.position="none", axis.title.y=element_blank(), plot.title=element_text(hjust=.5, face="bold")) +
    xlab("Hazard Ratio") + ggtitle(paste0(gene, ttl))
}

## ---- 5. main loop: set x cohort -------------------------------------
datasets <- list(Males=tcga_m, Females=tcga_f, All=tcga_all)
for (set in names(genesets)) {
  genes <- genesets[[set]]
  sub <- file.path(outtop, set); pf <- file.path(sub,"forest_plots"); kf <- file.path(sub,"km_plots")
  for (d in c(sub,pf,kf)) dir.create(d, showWarnings=FALSE, recursive=TRUE)
  wide_list <- list()
  for (coh in names(datasets)) {
    ph <- datasets[[coh]]
    suffix <- paste0(substr(set,1,1), substr(coh,1,1))   # MM/MF/MA, FM/FF/FA
    lng <- cox_long(ph, genes)
    write.csv(lng, file.path(sub, paste0("cox_long_",set,"_",coh,".csv")), row.names=FALSE)
    wd  <- cox_wide(ph, genes, suffix)
    write.csv(wd,  file.path(sub, paste0("cox_wide_",set,"_",coh,".csv")), row.names=FALSE)
    wide_list[[coh]] <- if (nrow(wd)) dplyr::rename(wd, !!paste0("Overall_pvalue_",coh) := Overall_pvalue) else NULL
    ## forest plots
    for (g in genes) tryCatch({
      p <- forest_plot(g, lng, paste0(" | ",set," | ",coh))
      if (!is.null(p)) ggsave(file.path(pf, paste0("forest_",g,"_",coh,".pdf")), p, width=6, height=4)
    }, error=function(e) NULL)
    ## KM plots (median split)
    if (have_surv) for (g in genes) {
      bc <- paste0(g,"_binary"); if (!bc %in% colnames(ph)) next
      tryCatch({
        fit <- survminer::surv_fit(as.formula(paste0("Surv(",OS_var,",",OS_event,") ~ ",bc)), data=ph)
        gp <- survminer::ggsurvplot(fit, data=ph, pval=TRUE, risk.table=TRUE, legend.title=g,
                                    legend.labs=c("Low","High"), xlab="OS (months)", palette="npg",
                                    ggtheme=theme_bw(), title=paste0(g," - ",set," | ",coh))
        ggsave(file.path(kf, paste0("KM_",g,"_",coh,".pdf")), gp$plot, width=7, height=5)
      }, error=function(e) NULL)
    }
    cat("  ", set, "x", coh, ":", length(genes), "TFs\n")
  }
  ## COMBINED across cohorts
  wide_list <- Filter(Negate(is.null), wide_list)
  combined <- Reduce(function(a,b) full_join(a,b, by="SYMBOL_UPPER"), wide_list)
  write.csv(combined, file.path(sub, paste0("cox_wide_",set,"_COMBINED.csv")), row.names=FALSE)
  ## pretty allTFs table: Motif + Expression HR/CI/SE/p (male & female patients)
  si <- substr(set,1,1)
  pretty <- combined %>% transmute(
    Motif = tf2motif(set)$motif.name[match(SYMBOL_UPPER, tf2motif(set)$TF_up)],
    TF = SYMBOL_UPPER,
    `Multivariate p value male`            = .data[[paste0("Overall_pvalue_Males")]],
    `Gene expression HR male`              = .data[[paste0("Expression_coef_",si,"M")]],
    `Gene expression HR lower CI male`     = .data[[paste0("Expression_Lower95_",si,"M")]],
    `Gene expression HR upper CI male`     = .data[[paste0("Expression_Upper95_",si,"M")]],
    `Gene expression HR standard error male` = .data[[paste0("Expression_SE_",si,"M")]],
    `Gene expression p value male`         = .data[[paste0("Expression_pvalue_",si,"M")]],
    `Multivariate p value female`          = .data[[paste0("Overall_pvalue_Females")]],
    `Gene expression HR female`            = .data[[paste0("Expression_coef_",si,"F")]],
    `Gene expression HR lower CI female`   = .data[[paste0("Expression_Lower95_",si,"F")]],
    `Gene expression HR upper CI female`   = .data[[paste0("Expression_Upper95_",si,"F")]],
    `Gene expression HR standard error female` = .data[[paste0("Expression_SE_",si,"F")]],
    `Gene expression p value female`       = .data[[paste0("Expression_pvalue_",si,"F")]])
  if (have_xlsx) openxlsx::write.xlsx(pretty, file.path(sub, paste0("cox_wide_",set,"_COMBINED_allTFs.xlsx")))
  write.csv(pretty, file.path(sub, paste0("cox_wide_",set,"_COMBINED_allTFs.csv")), row.names=FALSE)
  cat(set, "done ->", sub, "\n")
}
cat("\nAll done. Output root:\n  ", outtop, "\n")
