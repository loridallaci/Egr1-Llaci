## =====================================================================
## Supplementary Fig. 2 — ALL panels rebuilt on the PAPER'S cells, using the
## authors' cell-type labels. Runs on .181 (seurat5 env) over the withMotifs_V4
## objects, subset to the paper's QC cells (barcode reconciliation) with sex +
## author_cell_type transferred from the CELLxGENE annotation.
##
## Produces (OUT):
##   cellnumbers_paperCells.pdf                  (a: nuclei by stage x sex)
##   RNA_UMAP_paperCells_bySex.pdf / _byCelltype.pdf   (b: RPCA, integrated M/F)
##   ATAC_UMAP_paperCells_bySex.pdf / _byCelltype.pdf  (b: rLSI, integrated M/F)
##   coords_<stage>.csv                          (barcode, sex, celltype, umap coords)
##   chromvar_input_<stage>.rds                  (peaks counts + ranges + EGR1 col + meta)
##                                                 -> transferred to Windows for chromVAR (b2)
##
## Run:  LD_LIBRARY_PATH unneeded; use the seurat5 env with forced libpath:
##   ~/miniconda3/envs/seurat5/bin/Rscript -e '.libPaths(".../seurat5/lib/R/library");
##     source("Figure_SupplFig2_paperCells_panels.R")'
## =====================================================================
suppressMessages({ library(Seurat); library(Signac); library(ggplot2); library(patchwork); library(ggrepel); library(GenomicRanges) })
options(future.globals.maxSize = 400 * 1024^3); set.seed(1)

data_dir  <- "/home/lllaci/data/cortex_development"
out       <- "/home/lllaci/cortex_paperCells_panels"
annot_csv <- "/home/lllaci/paper_cell_annotation.csv"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

stage_levels <- c("LateFetal","Infant","Child","Adol","Adult")
stage_file <- c(LateFetal="LaFet", Infant="Inf", Child="Child", Adol="Adol", Adult="Adult")
stage_donor <- list(LateFetal=c("1"="LaFet1","2"="LaFet2"), Infant=c("1"="Inf1","2"="Inf2"),
                    Child=c("1"="Child1","2"="Child2"), Adol=c("1"="Adol1","2"="Adol2"),
                    Adult=c("1"="Adult1","2"="Adult2"))
sex_cols <- c(male="#4A6FE3", female="#F39AC9")
ct_order <- c("RG","IPC","EN-fetal-early","EN-fetal-late","EN","IN-fetal","IN-MGE","IN-CGE",
              "OPC","Oligodendrocytes","Astrocytes","Microglia","Endothelial","Pericytes","VSMC")
ct_cols <- c(RG="#8C564B", IPC="#E377C2", "EN-fetal-early"="#AEC7E8", "EN-fetal-late"="#3B6DB3",
             EN="#1F77B4", "IN-fetal"="#FF9896", "IN-MGE"="#D62728", "IN-CGE"="#E45756",
             OPC="#9467BD", Oligodendrocytes="#2CA02C", Astrocytes="#FF7F0E", Microglia="#17BECF",
             Endothelial="#7F7F7F", Pericytes="#4D4D4D", VSMC="#C49C94")
big_theme <- theme_classic(base_size = 16) +
  theme(plot.title=element_text(size=18,face="bold"), axis.title=element_text(size=16),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.title=element_text(size=14), legend.text=element_text(size=13))

paper <- read.csv(annot_csv, stringsAsFactors=FALSE); paper$key <- paste(paper$donor, paper$core, sep=":")

subset_to_paper <- function(obj, stage) {
  bc <- colnames(obj); core <- sub("-[0-9]+$","",bc); suf <- sub(".*-","",bc)
  key <- paste(stage_donor[[stage]][suf], core, sep=":"); m <- match(key, paper$key); keep <- !is.na(m)
  cat(sprintf("  %s: local %d -> paper-matched %d (%.0f%%)\n", stage, length(bc), sum(keep), 100*mean(keep)))
  obj <- obj[, keep]; obj$sex <- paper$sex[m[keep]]
  obj$celltype <- factor(paper$celltype[m[keep]], levels = ct_order); obj
}

coords_all <- list(); counts_list <- list()
for (st in stage_levels) {
  cat("\n#### ", st, " ####\n")
  fp <- Sys.glob(file.path(data_dir, paste0(stage_file[[st]], "_*withMotifs_V4_122825.rds")))[1]
  o <- readRDS(fp)
  for (a in intersect(c("peaks","ATAC"), Assays(o))) suppressWarnings(try(Fragments(o[[a]]) <- NULL, silent=TRUE))
  o <- subset_to_paper(o, st)

  ## ---- export chromVAR inputs (peaks counts + ranges + EGR1 motif col) ----
  pk <- o[["peaks"]]
  md <- Motifs(pk); mm <- md@data                       # peaks x motifs
  egr1_id <- names(which(unlist(md@motif.names) == "EGR1"))[1]
  if (is.na(egr1_id)) egr1_id <- "MA0162.4"
  egr1_col <- as.numeric(mm[, egr1_id])
  saveRDS(list(counts = GetAssayData(pk, layer="counts"),
               ranges = as.data.frame(granges(pk))[,c("seqnames","start","end")],
               egr1 = egr1_col, egr1_id = egr1_id,
               meta = data.frame(barcode=colnames(o), sex=o$sex, celltype=as.character(o$celltype))),
          file.path(out, paste0("chromvar_input_", st, ".rds")))

  ## ---- RNA UMAP (RPCA, integrated across sex) ----
  DefaultAssay(o) <- "RNA"
  if (!inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- as(o[["RNA"]], "Assay5")
  o[["RNA"]] <- split(o[["RNA"]], f = o$sex)
  o <- NormalizeData(o, verbose=FALSE); o <- FindVariableFeatures(o, verbose=FALSE)
  o <- ScaleData(o, verbose=FALSE); o <- RunPCA(o, verbose=FALSE)
  o <- IntegrateLayers(o, method=RPCAIntegration, orig.reduction="pca", new.reduction="integrated.RPCA", verbose=FALSE)
  o <- RunUMAP(o, reduction="integrated.RPCA", dims=1:30, reduction.name="umap.rna", verbose=FALSE)
  o[["RNA"]] <- JoinLayers(o[["RNA"]])
  er <- Embeddings(o, "umap.rna")

  ## ---- ATAC UMAP (rLSI, integrated across sex) ----
  DefaultAssay(o) <- "peaks"
  o <- RunTFIDF(o, verbose=FALSE); o <- FindTopFeatures(o, min.cutoff=5); o <- RunSVD(o, verbose=FALSE)
  ol <- SplitObject(o, split.by="sex")
  ol <- lapply(ol, function(x){ x<-RunTFIDF(x,verbose=FALSE); x<-FindTopFeatures(x,min.cutoff=5); x<-RunSVD(x,verbose=FALSE); x })
  anch <- FindIntegrationAnchors(object.list=ol, anchor.features=rownames(o), reduction="rlsi", dims=2:30, verbose=FALSE)
  integ <- IntegrateEmbeddings(anchorset=anch, reductions=o[["lsi"]], new.reduction.name="integrated_lsi", dims.to.integrate=1:30, verbose=FALSE)
  integ <- RunUMAP(integ, reduction="integrated_lsi", dims=2:30, reduction.name="umap.atac", verbose=FALSE)
  ea <- Embeddings(integ, "umap.atac")

  df <- data.frame(barcode=colnames(o), sex=factor(o$sex,levels=c("male","female")),
                   celltype=factor(o$celltype, levels=ct_order),
                   rna_x=er[colnames(o),1], rna_y=er[colnames(o),2],
                   atac_x=ea[colnames(o),1], atac_y=ea[colnames(o),2])
  write.csv(df, file.path(out, paste0("coords_", st, ".csv")), row.names=FALSE)
  coords_all[[st]] <- df; counts_list[[st]] <- as.data.frame(table(sex=df$sex)); counts_list[[st]]$stage <- st
  rm(o, integ, ol, anch); gc()
}

## ================= PANELS =================
present_all <- ct_order[ct_order %in% unlist(lapply(coords_all, function(d) as.character(d$celltype)))]
ct_scale <- function() scale_colour_manual(values=ct_cols, limits=present_all, drop=FALSE, name="Cell type (paper)",
                        guide=guide_legend(override.aes=list(size=3,alpha=1), ncol=1))

mk_sex <- function(d, ttl, xx, yy) ggplot(d[sample(nrow(d)),], aes(.data[[xx]], .data[[yy]], colour=sex)) +
  geom_point(size=0.35, alpha=0.7) + scale_colour_manual(values=sex_cols, name="Sex",
    guide=guide_legend(override.aes=list(size=3,alpha=1))) + coord_fixed() + labs(title=ttl,x=NULL,y=NULL) + big_theme
mk_ct <- function(d, ttl, xx, yy) { cent <- aggregate(cbind(get(xx),get(yy))~celltype, d, median)
  names(cent) <- c("celltype","cx","cy")
  ggplot(d[sample(nrow(d)),], aes(.data[[xx]], .data[[yy]], colour=celltype)) + geom_point(size=0.35, alpha=0.75) +
  geom_text_repel(data=cent, aes(cx,cy,label=celltype), inherit.aes=FALSE, size=4, fontface="bold", seed=1,
    box.padding=0.4, min.segment.length=0, colour="black", bg.color="white", bg.r=0.15, max.overlaps=Inf) +
  coord_fixed() + labs(title=ttl,x=NULL,y=NULL) + big_theme }

save_row <- function(plots, file, w=5.2*5+2.5) {
  m <- wrap_plots(plots, nrow=1)
  ggsave(file.path(out, paste0(file,".pdf")), m, width=w, height=6.0, limitsize=FALSE)
  ggsave(file.path(out, paste0(file,".png")), m, width=w, height=6.0, dpi=150, limitsize=FALSE) }

# legend only on last panel to keep one legend
withleg <- function(fn, colscale) lapply(seq_along(stage_levels), function(i){
  d <- coords_all[[stage_levels[i]]]; p <- fn(d, stage_levels[i])
  p + colscale + theme(legend.position = if (i==length(stage_levels)) "right" else "none") })

save_row(withleg(function(d,t) mk_sex(d,t,"rna_x","rna_y"), scale_colour_manual(values=sex_cols,name="Sex",guide=guide_legend(override.aes=list(size=3,alpha=1)))), "RNA_UMAP_paperCells_bySex")
save_row(withleg(function(d,t) mk_ct (d,t,"rna_x","rna_y"), ct_scale()), "RNA_UMAP_paperCells_byCelltype")
save_row(withleg(function(d,t) mk_sex(d,t,"atac_x","atac_y"), scale_colour_manual(values=sex_cols,name="Sex",guide=guide_legend(override.aes=list(size=3,alpha=1)))), "ATAC_UMAP_paperCells_bySex")
save_row(withleg(function(d,t) mk_ct (d,t,"atac_x","atac_y"), ct_scale()), "ATAC_UMAP_paperCells_byCelltype")

## ---- cell-number bar ----
cdf <- do.call(rbind, counts_list); cdf$stage <- factor(cdf$stage, levels=stage_levels)
write.csv(cdf[,c("stage","sex","Freq")], file.path(out,"CellCounts_paperCells_byStageSex.csv"), row.names=FALSE)
pc <- ggplot(cdf, aes(stage, Freq, fill=sex)) + geom_col(position=position_dodge(0.75), width=0.7, colour="black") +
  scale_fill_manual(values=sex_cols, name="Sex") +
  labs(title="Cortex development (paper's cells): nuclei by stage and sex", x="developmental stage", y="number of nuclei") +
  theme_classic(base_size=16) + theme(plot.title=element_text(size=18,face="bold"), axis.title=element_text(size=16),
    axis.text=element_text(size=14), axis.text.x=element_text(angle=15,hjust=1), legend.position="top",
    legend.title=element_text(size=14), legend.text=element_text(size=14))
ggsave(file.path(out,"cellnumbers_paperCells.pdf"), pc, width=8, height=5.5)
ggsave(file.path(out,"cellnumbers_paperCells.png"), pc, width=8, height=5.5, dpi=150)
cat("\nDONE ->", out, "\n"); print(cdf[,c("stage","sex","Freq")])
