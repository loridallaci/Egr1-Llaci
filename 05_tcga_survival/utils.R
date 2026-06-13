# =============================================================================
# 05_tcga_survival: Utility Functions
# =============================================================================
# Description:
#   Helper functions used across all TCGA survival analysis scripts.
#   Source this file at the top of each script:
#   source("05_tcga_survival/utils.R")
#
# Functions:
#   - Coxph.HR.CI()          : Extract HR and 95% CI from coxph object
#   - Format.Prob.CI()       : Format HR and CI as readable string
#   - smed()                 : Extract median survival from survfit object
#   - MakeSurvPlot()         : Generate Kaplan-Meier plot with Cox annotation
#   - run_cox_for_group()    : Run multivariate Cox regression for a gene list
#   - create_forest_plot()   : Generate publication-ready forest plot
#   - run_permutation_test() : Permutation test for independent prognostic value
# =============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# --- Cox HR and CI extraction -------------------------------------------------

Coxph.HR.CI <- function(coxph.obj) {
  
  check.assump <- cox.zph(coxph.obj, transform = "identity")
  ProportionalAssumption.p <- check.assump$table[, "p"]
  
  if (class(coxph.obj)[1] == "coxph") {
    res0     <- summary(coxph.obj)
    coef.mat <- res0$coef
    hr       <- coef.mat[, "exp(coef)"]
    pval     <- coef.mat[, "Pr(>|z|)"]
    se.coef  <- if ("robust se" %in% colnames(coef.mat)) {
      coef.mat[, "robust se"]
    } else {
      coef.mat[, "se(coef)"]
    }
    concord <- res0$concordance
    
  } else {
    stop("Error: need a coxph object!")
  }
  
  low.CI   <- hr * exp(-1.96 * se.coef)
  upper.CI <- hr * exp(1.96  * se.coef)
  
  out0 <- data.frame(
    HR             = round(hr, 4),
    lower          = round(low.CI, 4),
    upper          = round(upper.CI, 4),
    Pvalue         = pval,
    SE.coef        = se.coef,
    concordance    = concord[1],
    concordance.SE = concord[2]
  )
  rownames(out0) <- rownames(coef.mat)
  out0$"95%CI"  <- Format.Prob.CI(
    data.frame(HR = hr, lower = low.CI, upper = upper.CI),
    dec = 2, sep.char = "~"
  )
  out0 <- data.frame(Level = rownames(out0), out0,
                     stringsAsFactors = FALSE, check.names = FALSE)
  return(out0)
}

# --- Format CI as string ------------------------------------------------------

Format.Prob.CI <- function(dat0, dec = 2, sep.char = "~") {
  out <- apply(dat0, 1, function(zz) {
    zz   <- round(as.numeric(zz), dec)
    temp <- paste0(zz[1], " (", paste(zz[2], zz[3], sep = sep.char), ")")
    as.matrix(temp)
  })
  as.matrix(out)
}

# --- Median survival extractor ------------------------------------------------

smed <- function(x) {
  dat0 <- summary(x)$table
  if (!is.null(rownames(dat0))) {
    data.frame(Variable = rownames(dat0), dat0, check.names = FALSE)
  } else {
    dat0
  }
}

# --- Kaplan-Meier plot with Cox annotation ------------------------------------

MakeSurvPlot <- function(PFS = "PFS", event = "PFS_event", geneID = "Gene",
                         dat, save = FALSE,
                         legend.labs = c("Low", "High"),
                         legend.title = "gene", ...) {
  
  require("survival")
  require("survminer")
  
  dat  <- na.exclude(dat[, c(PFS, event, geneID)])
  fit  <- surv_fit(
    as.formula(paste("Surv(", PFS, ",", event, ") ~", geneID)),
    data = dat
  )
  
  med0 <- smed(fit)
  if (is.null(legend.labs)) legend.labs <- levels(factor(dat[, geneID]))
  med0$Variable <- legend.labs
  
  legend0 <- paste0(
    med0$Variable, ": n=", med0$n.start,
    " event=", med0$events,
    " med=", round(med0$median, 2),
    " (", round(med0$"0.95LCL", 2), "~", round(med0$"0.95UCL", 2), ")"
  )
  
  coxfit0 <- coxph(
    as.formula(paste("Surv(", PFS, ",", event, ") ~", geneID)),
    data = dat
  )
  HR <- Coxph.HR.CI(coxfit0)$"95%CI"
  
  plot1 <- ggsurvplot(
    fit,
    data             = dat,
    tables.theme     = theme_cleantable(),
    risk.table       = TRUE,
    pval             = TRUE,
    conf.int         = FALSE,
    xlab             = "OS (months)",
    legend.title     = legend.title,
    legend.labs      = legend.labs,
    legend           = "top",
    surv.median.line = "hv",
    ggtheme          = theme_bw(),
    risk.table.y.text.col = TRUE,
    risk.table.y.text     = TRUE,
    font.x           = c(12, "bold"),
    font.y           = c(12, "bold"),
    font.tickslab    = c(10, "bold"),
    font.main        = c(12, "bold"),
    font.legend      = c(15, "bold"),
    pval.size        = 4,
    pval.coord       = c(30, 0.55),
    palette          = "npg",
    ...
  )
  
  plot1$plot <- plot1$plot +
    ggplot2::annotate("text", x = 30, y = 0.99, hjust = 0,
                      label = legend0[1], color = "#FF2700") +
    ggplot2::annotate("text", x = 30, y = 0.85, hjust = 0,
                      label = legend0[2], color = "#008FD5") +
    ggplot2::annotate("text", x = 30, y = 0.7,  hjust = 0,
                      label = paste0("HR (high vs. low) = ", HR))
  
  print(plot1)
  return(plot1$plot)
}

# --- Run multivariate Cox regression for a gene list -------------------------
# Returns a tidy data frame with HR, CI, SE, and p-values per gene/variable

run_cox_for_group <- function(pheno_data, genes,
                              OS_var     = "survival",
                              OS_event   = "status",
                              covariates = c("Recurrence", "Age",
                                             "Subtype", "MGMT_status")) {
  
  survival_table <- data.frame()
  
  for (gene_name in genes) {
    
    if (!gene_name %in% colnames(pheno_data)) {
      warning(paste("Gene", gene_name, "not found in data. Skipping."))
      next
    }
    
    formula <- as.formula(
      paste("Surv(", OS_var, ",", OS_event, ") ~", gene_name, "+",
            paste(covariates, collapse = "+"))
    )
    
    tryCatch({
      cox_obj <- coxph(formula, data = pheno_data)
      cox_sum <- summary(cox_obj)
      
      df <- as.data.frame(cox_sum$coefficients)
      ci <- as.data.frame(cox_sum$conf.int)
      
      result <- data.frame(
        Gene     = gene_name,
        Variable = rownames(df),
        HR       = ci$`exp(coef)`,
        Lower95  = ci$`lower .95`,
        Upper95  = ci$`upper .95`,
        SE       = df$`se(coef)`,
        Pvalue   = df$`Pr(>|z|)`
      )
      
      survival_table <- rbind(survival_table, result)
      
    }, error = function(e) {
      warning(paste("Error with gene", gene_name, ":", e$message))
    })
  }
  
  # Standardize variable name labels
  survival_table$Variable <- gsub("MGMT_statusUnmethylated",
                                  "MGMT (unmethylated)",    survival_table$Variable)
  survival_table$Variable <- gsub("RecurrenceRecurrent",
                                  "Recurrence (recurrent)", survival_table$Variable)
  survival_table$Variable <- gsub("RecurrenceSecondary",
                                  "Recurrence (secondary)", survival_table$Variable)
  survival_table$Variable <- gsub("SubtypeMesenchymal",
                                  "Mesenchymal Subtype",    survival_table$Variable)
  survival_table$Variable <- gsub("SubtypeProneural",
                                  "Proneural Subtype",      survival_table$Variable)
  
  return(survival_table)
}

# --- Forest plot for a single gene --------------------------------------------

create_forest_plot <- function(gene_name, survival_data,
                               title_suffix = "",
                               var_order    = NULL,
                               xlim_range   = c(0, 2)) {
  
  gene_data   <- survival_data %>% dplyr::filter(Gene == gene_name)
  unique_vars <- unique(gene_data$Variable)
  
  if (is.null(var_order)) {
    var_order <- c(gene_name,
                   "Age",
                   "MGMT (unmethylated)",
                   "Recurrence (recurrent)",
                   "Recurrence (secondary)",
                   "Mesenchymal Subtype",
                   "Proneural Subtype")
  }
  var_order <- var_order[var_order %in% unique_vars]
  
  plot_df <- gene_data %>%
    mutate(
      signif   = ifelse(Pvalue <= 0.05, "p<=0.05", "ns"),
      Variable = factor(Variable, levels = rev(var_order))
    )
  
  ggplot(plot_df, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = signif), size = 4) +
    geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.3) +
    scale_color_manual(values = c("p<=0.05" = "red", "ns" = "black")) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.y    = element_blank(),
      axis.text.y     = element_text(size = 12),
      plot.title      = element_text(hjust = 0.5, face = "bold")
    ) +
    xlab("Hazard Ratio") +
    ggtitle(paste0(gene_name, title_suffix)) +
    xlim(xlim_range)
}

# --- Permutation test for independent prognostic value ------------------------
# For each gene, permutes expression n_perm times and compares the observed
# log-HR to the null distribution. A significant perm_pvalue means the gene
# adds prognostic value INDEPENDENT of the covariates.
#
# Arguments:
#   pheno_data  : phenotype data frame (e.g. tcga_pheno_male)
#   genes       : character vector of gene names to test
#   obs_results : output of run_cox_for_group() for the same pheno_data
#   OS_var      : survival time column name
#   OS_event    : event indicator column name
#   covariates  : covariate names (must match run_cox_for_group)
#   n_perm      : number of permutations (default 1000; use 5000 for publication)
#
# Returns a data frame with one row per gene containing:
#   Gene, obs_logHR, obs_pvalue, perm_pvalue, perm_padj (BH), perm_sig

run_permutation_test <- function(pheno_data,
                                 genes,
                                 obs_results,
                                 OS_var     = "survival",
                                 OS_event   = "status",
                                 covariates = c("Recurrence", "Age",
                                                "Subtype", "MGMT_status"),
                                 n_perm     = 1000) {
  
  results <- data.frame()
  
  for (gene_name in genes) {
    
    if (!gene_name %in% colnames(pheno_data)) {
      warning(paste("Gene", gene_name, "not found in data. Skipping."))
      next
    }
    
    # Observed log-HR for this gene from the real Cox model
    obs_row <- obs_results %>%
      filter(Gene == gene_name, Variable == gene_name)
    
    if (nrow(obs_row) == 0) {
      warning(paste("No observed result for gene", gene_name, "- skipping."))
      next
    }
    
    obs_logHR  <- log(obs_row$HR[1])
    obs_pvalue <- obs_row$Pvalue[1]
    
    formula <- as.formula(
      paste("Surv(", OS_var, ",", OS_event, ") ~", gene_name, "+",
            paste(covariates, collapse = "+"))
    )
    
    # Permutation null distribution
    perm_logHRs <- numeric(n_perm)
    
    for (p in seq_len(n_perm)) {
      perm_data <- pheno_data
      perm_data[[gene_name]] <- sample(pheno_data[[gene_name]])
      
      tryCatch({
        perm_cox <- coxph(formula, data = perm_data)
        perm_logHRs[p] <- perm_cox$coefficients[gene_name]
      }, error = function(e) {
        perm_logHRs[p] <<- NA_real_
      })
    }
    
    perm_logHRs <- perm_logHRs[!is.na(perm_logHRs)]
    
    # Two-sided permutation p-value
    perm_pvalue <- mean(abs(perm_logHRs) >= abs(obs_logHR))
    
    results <- rbind(results, data.frame(
      Gene        = gene_name,
      obs_logHR   = obs_logHR,
      obs_HR      = obs_row$HR[1],
      obs_pvalue  = obs_pvalue,
      perm_pvalue = perm_pvalue,
      n_perm_used = length(perm_logHRs),
      stringsAsFactors = FALSE
    ))
  }
  
  # BH correction across all genes
  results$perm_padj <- p.adjust(results$perm_pvalue, method = "BH")
  results$perm_sig  <- ifelse(results$perm_padj <= 0.05, "Yes", "No")
  
  return(results)
}
