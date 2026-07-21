# ============================================================================
# Cortex development — rebuild per-stage objects on the PAPER's cells (100%)
# ============================================================================
# Supersedes the QC-filtered build (cortex_dev_build_objects.R) which only
# recovered ~90% of the authors' cells. Here we start from the RAW ARC matrix
# (raw_feature_bc_matrix.h5) and keep EXACTLY the paper's barcodes — validated
# to recover 100.0% of the paper's cells (30,769/30,769 across LaFet/Inf/
# Child/Adult). Peaks + annotation are REUSED from the existing withMotifs_V4
# objects (same peak set as the current figure; only the cell set grows), so no
# MACS2/EnsDb needed. Output feeds cortex_dev_add_motifs.R unchanged.
#
# Run on .181 in the RENIN env (Seurat 4.1.1 / Signac 1.10, matches the objects):
#   LD_LIBRARY_PATH=$HOME/local/lib /home/lllaci/R-4.2.2/bin/Rscript cortex_dev_build_paperCells.R
# (Adol is added later, once its cellranger-arc run finishes.)
# ============================================================================
suppressMessages({library(Seurat); library(Signac); library(Matrix)})
set.seed(1234)

arc_dir  <- "/home/lllaci/data/cortex_development_arc"        # raw/filtered h5 + fragments per donor
obj_dir  <- "/home/lllaci/data/cortex_development"            # existing withMotifs_V4 objects (peaks source)
out_dir  <- "/home/lllaci/data/cortex_development_paperCells"; dir.create(out_dir, showWarnings=FALSE)
paper    <- read.csv("/home/lllaci/paper_cell_annotation.csv", stringsAsFactors=FALSE)

# stage -> 2 donors, in the order that fixes suffix -1 = donor1, -2 = donor2
stages <- list(LateFetal=c("LaFet1","LaFet2"), Infant=c("Inf1","Inf2"),
               Child=c("Child1","Child2"),     Adult=c("Adult1","Adult2"))
src_obj <- c(LateFetal="LaFet", Infant="Inf", Child="Child", Adult="Adult")   # existing-object prefixes
wm <- function(p) file.path(obj_dir, paste0(p,"_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable_withMotifs_V4_122825.rds"))

# ---- build one donor: RNA (raw h5, paper cells) + a fragment object with cell map ----
build_donor <- function(donor, suffix) {
  cts  <- Read10X_h5(file.path(arc_dir, donor, "raw_feature_bc_matrix.h5"))
  rna  <- cts[["Gene Expression"]]
  core <- sub("-1$", "", colnames(rna))
  keep <- core %in% paper$core[paper$donor == donor]
  rna  <- rna[, keep]
  old  <- colnames(rna)                          # fragment-file barcodes ({core}-1)
  new  <- paste0(core[keep], "-", suffix)        # object barcodes ({core}-1 or -2)
  colnames(rna) <- new
  # Signac cells map: names = object barcode, value = fragment-file barcode
  frag <- CreateFragmentObject(path = file.path(arc_dir, donor, "atac_fragments.tsv.gz"),
                               cells = setNames(old, new))
  list(rna = rna, frag = frag, n = ncol(rna))
}

build_stage <- function(stage) {
  cat("\n=================", stage, "=================\n")
  donors <- stages[[stage]]
  d1 <- build_donor(donors[1], 1); d2 <- build_donor(donors[2], 2)
  cat(sprintf("  %s %d + %s %d = %d paper cells\n", donors[1], d1$n, donors[2], d2$n, d1$n + d2$n))

  # --- RNA object (merge the two donors) ---
  rna <- cbind(d1$rna, d2$rna)
  obj <- CreateSeuratObject(counts = rna, assay = "RNA")

  # --- reuse the existing stage peak set + annotation ---
  ex     <- readRDS(wm(src_obj[[stage]]))
  peaks  <- granges(ex[["peaks"]])
  annot  <- Annotation(ex[["peaks"]])
  rm(ex); gc()

  # --- quantify those peaks on the paper cells over both donors' fragments ---
  fmat <- FeatureMatrix(fragments = list(d1$frag, d2$frag),
                        features = peaks, cells = colnames(obj))
  obj[["peaks"]] <- CreateChromatinAssay(counts = fmat,
                        fragments = list(d1$frag, d2$frag), annotation = annot)

  # --- sex + celltype straight from the paper annotation (donor:core key) ---
  bc    <- colnames(obj); core <- sub("-[12]$", "", bc)
  donor <- ifelse(grepl("-1$", bc), donors[1], donors[2])
  m     <- match(paste(donor, core, sep=":"), paste(paper$donor, paper$core, sep=":"))
  obj$sex      <- paper$sex[m]
  obj$celltype <- paper$celltype[m]
  obj$stage    <- stage
  cat(sprintf("  sex: %s | NA sex: %d\n",
              paste(names(table(obj$sex)), table(obj$sex), collapse=" "), sum(is.na(obj$sex))))

  out <- file.path(out_dir, paste0(src_obj[[stage]],
           "_aggregated_object_Roussos_112223_RNAandPeakAssaysAvailable.rds"))
  saveRDS(obj, out); cat("  saved ->", out, "\n")
}

for (st in names(stages)) build_stage(st)
cat("\n=== Done. Paper-cell objects in", out_dir, "===\n")
cat("Next: point cortex_dev_add_motifs.R at", out_dir, "-> withMotifs_V4 -> Supp Fig 2 panels.\n")
