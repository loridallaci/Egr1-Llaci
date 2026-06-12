suppressMessages(library(Seurat)) 
suppressMessages(library(Signac))
suppressMessages(library(SeuratWrappers)) 
suppressMessages(library(RENIN)) 
suppressMessages(library(harmony))
#plan(sequential)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)

# Work on ChromVar object
# original (author's machine): '/home/lllaci/data/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds'
lot6 <- readRDS('output/female_male_aggregated_081722_seuratObject_multiome_081721_filtered_012925_withPeaks_chromVARadded_111425.rds')
obj <- lot6
DefaultAssay(obj) <- 'RNA'
# Normalize RNA data
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Verify the range changed
range(GetAssayData(obj, assay = "RNA", slot = "data")["Cd44", ])
range(GetAssayData(obj, assay = "RNA", slot = "data")["Crabp1", ])

Idents(obj) <- "sex"

DefaultAssay(obj) <- 'ATAC'

# remove old path to fragments
Fragments(obj[["ATAC"]]) <- NULL
Fragments(obj[["ATAC"]])
list()


# original (author's machine): "/home/lllaci/data/atac_fragments.tsv.gz"
fragpath <- "data/atac_fragments.tsv.gz"

frag.obj <- CreateFragmentObject(
  path = fragpath,
  cells = colnames(obj)   # ensures only your cells are used
)
Fragments(obj[["ATAC"]]) <- frag.obj


cairo_pdf("maleANDfemale_lot6_tracks_Crabp1_splitbysex.pdf",
          width = 6,
          height = 4)

p1 <- CoveragePlot(
  object = obj,
  region = "Crabp1",
  features = "Crabp1",
  expression.assay = "RNA",
  group.by = "sex",
  extend.upstream = 500,
  extend.downstream = 10000,
  cols = c("blue", "pink")
)

print(p1)

dev.off()


cairo_pdf("maleANDfemale_lot6_tracks_Igfbp2_splitbysex.pdf",
          width = 6,
          height = 4)

p1 <- CoveragePlot(
  object = obj,
  region = "Igfbp2",
  features = "Igfbp2",
  expression.assay = "RNA",
  group.by = "sex",
  extend.upstream = 500,
  extend.downstream = 10000,
  cols = c("blue", "pink")
)

print(p1)

dev.off()





