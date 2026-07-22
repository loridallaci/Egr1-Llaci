## Figure 5 (sterol story) — GSEA showing Egr1 KD coordinately de-represses the
## SREBP-2 / cholesterol program, more strongly in males than females.
## Independent MSigDB gene sets; ranks = DESeq2 Wald stat from the KD-vs-WT DE.
suppressMessages({library(fgsea); library(msigdbr); library(ggplot2)})
de_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_Egr1KD_gRNA3_vs_Neg1"
out    <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/09_drug_screen"
M  <- read.delim(file.path(de_dir,"Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt"))
Fm <- read.delim(file.path(de_dir,"Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"))

want <- c("HORTON_SREBF_TARGETS","HALLMARK_CHOLESTEROL_HOMEOSTASIS","REACTOME_CHOLESTEROL_BIOSYNTHESIS",
          "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP",
          "REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF",
          "WP_MEVALONATE_ARM_OF_CHOLESTEROL_BIOSYNTHESIS_PATHWAY",
          "GOBP_STEROL_BIOSYNTHETIC_PROCESS",
          "WP_STEROL_REGULATORY_ELEMENTBINDING_PROTEINS_SREBP_SIGNALING")
lab <- c(HORTON_SREBF_TARGETS="SREBF direct targets (Horton)",
         HALLMARK_CHOLESTEROL_HOMEOSTASIS="Cholesterol homeostasis (Hallmark)",
         REACTOME_CHOLESTEROL_BIOSYNTHESIS="Cholesterol biosynthesis (Reactome)",
         REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP="SREBP activates transcription (Reactome)",
         REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF="SREBP regulates cholesterol synth (Reactome)",
         WP_MEVALONATE_ARM_OF_CHOLESTEROL_BIOSYNTHESIS_PATHWAY="Mevalonate arm (WikiPathways)",
         GOBP_STEROL_BIOSYNTHETIC_PROCESS="Sterol biosynthetic process (GO:BP)",
         WP_STEROL_REGULATORY_ELEMENTBINDING_PROTEINS_SREBP_SIGNALING="SREBP signaling (WikiPathways)")
msig <- msigdbr(species="Homo sapiens")
sel  <- msig[msig$gs_name %in% want,]
pathways <- lapply(split(toupper(sel$gene_symbol), sel$gs_name), unique)

rankvec <- function(tab){ v<-tab$stat; names(v)<-toupper(tab$SYMBOL)
  v<-v[!is.na(v) & names(v)!="" & names(v)!="NA"]
  v<-tapply(v,names(v),function(x) x[which.max(abs(x))]); sort(v,decreasing=TRUE) }
rM<-rankvec(M); rF<-rankvec(Fm)
set.seed(1); fM<-fgsea(pathways,rM,minSize=8,maxSize=400,eps=0)
set.seed(1); fF<-fgsea(pathways,rF,minSize=8,maxSize=400,eps=0)

th <- theme_bw(base_size=16)+theme(
  axis.text=element_text(size=14,colour="black"), axis.title=element_text(size=16),
  plot.title=element_text(size=18,face="bold"), plot.subtitle=element_text(size=13,colour="grey30"),
  legend.text=element_text(size=14), legend.title=element_blank(), legend.position="top",
  plot.caption=element_text(size=12,hjust=0), panel.grid.minor=element_blank(),
  plot.title.position="plot")
star <- function(p) ifelse(is.na(p),"",ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",""))))

## ---- panel 1: male-vs-female NES bars ----
D <- rbind(data.frame(set=fM$pathway,sex="Male",  NES=fM$NES,padj=fM$padj),
           data.frame(set=fF$pathway,sex="Female",NES=fF$NES,padj=fF$padj))
ord <- fM$pathway[order(fM$NES)]
D$label <- factor(lab[D$set], levels=lab[ord])
D$sex   <- factor(D$sex, levels=c("Male","Female"))
D$txt   <- paste0(sprintf("%.1f",D$NES), star(D$padj))
p1 <- ggplot(D, aes(NES, label, fill=sex))+
  geom_col(position=position_dodge(width=.72), width=.66, colour="grey30")+
  geom_text(aes(label=txt), position=position_dodge(width=.72), hjust=-0.12, size=4.3)+
  scale_fill_manual(values=c(Male="#1874CD", Female="#C71585"), breaks=c("Male","Female"))+
  scale_x_continuous(limits=c(0,3.05), expand=c(0,0))+
  labs(title="Egr1 KD de-represses the SREBP-2\ncholesterol program (stronger in males)",
       subtitle="GSEA of independent MSigDB sets; NES > 0 = de-repressed after Egr1 KD",
       x="Normalized enrichment score (NES)", y=NULL,
       caption="* padj<0.05   ** <0.01   *** <0.001")+th
ggsave(file.path(out,"Figure5_SREBP_GSEA_NES.pdf"), p1, width=9.5, height=6.4, device=cairo_pdf)
ggsave(file.path(out,"Figure5_SREBP_GSEA_NES.png"), p1, width=9.5, height=6.4, dpi=300, bg="white")

## ---- panel 2: running-enrichment plot for the curated Horton SREBF targets (male) ----
i <- which(fM$pathway=="HORTON_SREBF_TARGETS")
ph <- plotEnrichment(pathways[["HORTON_SREBF_TARGETS"]], rM)+
  labs(title="SREBF direct targets (Horton) enrich\namong male Egr1-KD up-regulated genes",
       subtitle=sprintf("NES = %.2f,  padj = %.1e  (male KD)", fM$NES[i], fM$padj[i]),
       x="Gene rank in male KD  (up-regulated → down-regulated)", y="Running enrichment score")+th
ggsave(file.path(out,"Figure5_SREBP_GSEA_HORTON.pdf"), ph, width=7.5, height=5, device=cairo_pdf)
ggsave(file.path(out,"Figure5_SREBP_GSEA_HORTON.png"), ph, width=7.5, height=5, dpi=300, bg="white")

cat("wrote Figure5_SREBP_GSEA_NES.{pdf,png} and Figure5_SREBP_GSEA_HORTON.{pdf,png}\n")
