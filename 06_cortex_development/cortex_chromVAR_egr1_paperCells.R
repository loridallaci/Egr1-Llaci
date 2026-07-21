## =====================================================================
## Egr1 chromVAR motif activity on the PAPER-CELL (100%) cortex objects.
## Reads the chromvar_input_<stage>.rds exported by the panels script
## (counts, ranges, egr1 motif membership, meta), computes per-cell EGR1
## deviation z-scores, writes Egr1motif_activity_<stage>.csv, and draws the
## by-celltype x sex summary (panel b2 / extra panel).
## Runs on .181 in R-4.2.2 (Matrix < 1.6.4, so the chromVAR crossprod bug that
## forced the Windows old-Matrix recipe does not fire here). Defensive numeric
## coercions kept from Figures_SupplFig1e_chromVAR_Egr1_motif_fix.R just in case.
##   LD_LIBRARY_PATH=$HOME/local/lib /home/lllaci/R-4.2.2/bin/Rscript cortex_chromVAR_egr1_paperCells.R
## =====================================================================
suppressMessages({library(chromVAR); library(SummarizedExperiment); library(BSgenome.Hsapiens.UCSC.hg38)
                  library(Matrix); library(GenomicRanges); library(BiocParallel); library(ggplot2)})
register(SerialParam(), default = TRUE); set.seed(1)
setMethod("rowSums","CsparseMatrix",function(x,na.rm=FALSE,dims=1,...){if(!methods::is(x,"dMatrix"))x<-x*1;as.numeric(x%*%rep(1,ncol(x)))})
setMethod("colSums","CsparseMatrix",function(x,na.rm=FALSE,dims=1,...){if(!methods::is(x,"dMatrix"))x<-x*1;as.numeric(crossprod(rep(1,nrow(x)),x))})

dir <- "/home/lllaci/cortex_paperCells_panels_100pct"
stages <- c("LateFetal","Infant","Child","Adult")

for (st in stages) {
  cat("====", st, "====\n")
  inp <- readRDS(file.path(dir, paste0("chromvar_input_", st, ".rds")))
  counts <- as(as(inp$counts, "CsparseMatrix"), "dMatrix"); if (!is.double(counts@x)) counts@x <- as.double(counts@x)
  mm <- Matrix(as.numeric(inp$egr1 > 0), ncol = 1, sparse = TRUE); colnames(mm) <- "EGR1"
  mm <- as(as(mm, "CsparseMatrix"), "dMatrix"); if (!is.double(mm@x)) mm@x <- as.double(mm@x)
  gr <- makeGRangesFromDataFrame(inp$ranges)
  rse <- SummarizedExperiment(assays = list(counts = counts), rowRanges = gr)
  rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
  bg  <- getBackgroundPeaks(rse)
  dev <- computeDeviations(object = rse, annotations = mm, background_peaks = bg)
  z   <- as.numeric(deviationScores(dev)[1, ])
  out <- data.frame(barcode = inp$meta$barcode, sex = inp$meta$sex,
                    celltype = as.character(inp$meta$celltype), egr1_z = z, stage = st)
  write.csv(out, file.path(dir, paste0("Egr1motif_activity_", st, ".csv")), row.names = FALSE)
  cat(sprintf("  %s: %d cells, EGR1 z range %.2f..%.2f\n", st, nrow(out), min(z), max(z)))
}

## ---- pooled by-celltype x sex summary + plot ----
all <- do.call(rbind, lapply(stages, function(st) read.csv(file.path(dir, paste0("Egr1motif_activity_", st, ".csv")))))
summ <- do.call(rbind, lapply(unique(all$celltype), function(ct){
  s <- all[all$celltype==ct,]; m<-s$egr1_z[s$sex=="male"]; f<-s$egr1_z[s$sex=="female"]
  if(length(m)<3||length(f)<3) return(NULL)
  data.frame(celltype=ct, n_M=length(m), n_F=length(f), mean_M=mean(m), mean_F=mean(f),
             sem_M=sd(m)/sqrt(length(m)), sem_F=sd(f)/sqrt(length(f)),
             wilcox_p=suppressWarnings(wilcox.test(m,f)$p.value))}))
summ$p_adj <- p.adjust(summ$wilcox_p,"BH")
summ$star  <- cut(summ$p_adj,c(-Inf,1e-3,1e-2,5e-2,Inf),labels=c("***","**","*","ns"))
summ$overall <- (summ$mean_M*summ$n_M + summ$mean_F*summ$n_F)/(summ$n_M+summ$n_F)
summ <- summ[order(summ$overall),]; summ$celltype <- factor(summ$celltype, levels=summ$celltype)
write.csv(summ[order(-summ$overall),], file.path(dir,"EGR1_activity_byCelltypeSex_summary.csv"), row.names=FALSE)

long <- rbind(data.frame(celltype=summ$celltype,sex="male",z=summ$mean_M),
              data.frame(celltype=summ$celltype,sex="female",z=summ$mean_F))
long$sex <- factor(long$sex,levels=c("male","female"))
p <- ggplot(long, aes(z, celltype, colour=sex)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60") +
  geom_segment(data=summ, inherit.aes=FALSE, aes(x=mean_M,xend=mean_F,y=celltype,yend=celltype), colour="grey70", linewidth=0.8) +
  geom_point(size=4) +
  geom_text(data=data.frame(celltype=summ$celltype, star=summ$star, x=pmax(summ$mean_M,summ$mean_F)+0.06),
            inherit.aes=FALSE, aes(x=x,y=celltype,label=star), hjust=0, size=5, colour="grey25") +
  scale_colour_manual(values=c(male="#4A6FE3",female="#F39AC9"), name="Sex") +
  scale_x_continuous(expand=expansion(mult=c(0.03,0.12))) +
  labs(x="EGR1 motif activity (mean chromVAR z)", y=NULL,
       title="EGR1 motif activity by cell type and sex",
       subtitle="Developing human cortex, paper's cells (100%), 4 stages pooled") +
  theme_classic(base_size=16) +
  theme(plot.title=element_text(size=18,face="bold"), plot.subtitle=element_text(size=13),
        axis.title.x=element_text(size=16), axis.text=element_text(size=14,colour="black"),
        legend.text=element_text(size=14), legend.title=element_text(size=14), legend.position="top")
ggsave(file.path(dir,"SuppFig2_EGR1_activity_byCelltypeSex.pdf"), p, width=8, height=7)
ggsave(file.path(dir,"SuppFig2_EGR1_activity_byCelltypeSex.png"), p, width=8, height=7, dpi=300, bg="white")
cat("\nDone. per-cell CSVs + byCelltypeSex plot/summary in", dir, "\n")
