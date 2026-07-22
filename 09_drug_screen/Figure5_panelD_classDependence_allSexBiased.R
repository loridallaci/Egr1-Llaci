## Figure 5 panel D (alternative) — ALL sex-biased hits (both male- and female-biased),
## matching the panel C Venn: WT-only (27) | Shared (22) | KD-only (3).
## Same build as Figure5_panelD_classDependence.R but WITHOUT the male-biased (Q2<0) filter,
## so every q<=0.05 sex-difference hit is counted. NB: the 2 female-biased hits
## (Guanabenz, Prasterone) are both KD-only, so they add to the KD-only/"gained" segments.
suppressMessages({library(readxl); library(ggplot2)})
wb  <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Cell Titer Glo/Maxene Drug Screen/85compoundsCheck_analyses_July2026.xlsx"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/09_drug_screen"

d <- suppressMessages(read_excel(wb, sheet="Zc combined", col_names=FALSE))[-1,]
S <- data.frame(name=toupper(trimws(as.character(d[[9]]))),
                M_WT=as.numeric(d[[3]]), M_KD=as.numeric(d[[4]]),
                F_WT=as.numeric(d[[5]]), F_KD=as.numeric(d[[6]]))
q2 <- function(a,b)(a-b)/sqrt(2)
S$Q2_WT <- q2(S$M_WT,S$F_WT); S$Q2_KD <- q2(S$M_KD,S$F_KD)
S$sigWT <- p.adjust(2*pnorm(-abs(S$Q2_WT)),"BH")<=0.05
S$sigKD <- p.adjust(2*pnorm(-abs(S$Q2_KD)),"BH")<=0.05

cls <- c(
 "ATORVASTATIN CALCIUM"="Statin","MEVASTATIN"="Statin","LOVASTATIN"="Statin","SIMVASTATIN"="Statin",
 "PREDNISOLONE SODIUM PHOSPHATE"="Corticosteroid","FLUOCINOLONE ACETONIDE"="Corticosteroid","METHYLPREDNISOLONE"="Corticosteroid",
 "FLUMETHAZONE PIVALATE"="Corticosteroid","FLUDROCORTISONE ACETATE"="Corticosteroid","DEXAMETHASONE"="Corticosteroid",
 "BAZEDOXIFENE ACETATE"="Sex/steroid hormone","TOREMIFENE CITRATE"="Sex/steroid hormone","TAMOXIFEN CITRATE"="Sex/steroid hormone",
 "CLOMIPHENE CITRATE"="Sex/steroid hormone","ESTRONE"="Sex/steroid hormone","ESTRIOL"="Sex/steroid hormone",
 "ESTRADIOL CYPIONATE"="Sex/steroid hormone","ALLYLESTRENOL"="Sex/steroid hormone","NORETHINDRONE ACETATE"="Sex/steroid hormone",
 "PRASTERONE ACETATE"="Sex/steroid hormone","EPLERENONE"="Sex/steroid hormone",
 "ADEFOVIR DIPIVOXYL"="Antimetabolite","AZASERINE"="Antimetabolite","CARMOFUR"="Antimetabolite","FLOXURIDINE"="Antimetabolite",
 "AZACITIDINE"="Antimetabolite","MYCOPHENOLIC ACID"="Antimetabolite","THIOGUANINE"="Antimetabolite","FLUOROURACIL"="Antimetabolite","AZATHIOPRINE"="Antimetabolite",
 "AMSACRINE"="DNA-damaging","OXALIPLATIN"="DNA-damaging","MITOMYCIN"="DNA-damaging","AMINACRINE"="DNA-damaging",
 "DOCETAXEL"="Microtubule","VINBLASTINE SULFATE"="Microtubule",
 "FLUBENDAZOLE"="Antiparasitic","MONENSIN SODIUM"="Antiparasitic","NITARSONE"="Antiparasitic","ARTESUNATE"="Antiparasitic",
 "MEFLOQUINE HYDROCHLORIDE"="Antiparasitic","MEBENDAZOLE"="Antiparasitic","FENBENDAZOLE"="Antiparasitic","PYRVINIUM PAMOATE"="Antiparasitic","ATOVAQUONE"="Antiparasitic",
 "ITRACONAZOLE HYDROCHLORIDE"="Antimicrobial","BENZOXIQUINE"="Antimicrobial","OXICONAZOLE NITRATE"="Antimicrobial","SECNIDAZOLE"="Antimicrobial",
 "DIRITHROMYCIN"="Antimicrobial","OXYQUINOLINE HEMISULFATE"="Antimicrobial","FUSIDIC ACID"="Antimicrobial","FURAZOLIDONE"="Antimicrobial","BENZETHONIUM CHLORIDE"="Antimicrobial",
 "METHYLPHENIDATE HYDROCHLORIDE"="CNS/neuroactive","RILUZOLE"="CNS/neuroactive","METITEPINE MESYLATE"="CNS/neuroactive","KETANSERIN"="CNS/neuroactive",
 "PAROXETINE HYDROCHLORIDE"="CNS/neuroactive","SULOCTIDIL"="CNS/neuroactive","PIMOZIDE"="CNS/neuroactive","FLUPHENAZINE HYDROCHLORIDE"="CNS/neuroactive",
 "GUANABENZ ACETATE"="CNS/neuroactive","ETHOPROPAZINE HYDROCHLORIDE"="CNS/neuroactive",
 "BRONOPOL"="Oxidative/thiol","HYDROQUINONE"="Oxidative/thiol","MENADIONE"="Oxidative/thiol","PHENYLMERCURIC ACETATE"="Oxidative/thiol",
 "THIMEROSAL"="Oxidative/thiol","OLTIPRAZ"="Oxidative/thiol","ASCORBIC ACID"="Oxidative/thiol","CIANIDANOL [+-CATECHIN]"="Oxidative/thiol","ROTENONE"="Oxidative/thiol",
 "CELECOXIB"="NSAID","OXYPHENBUTAZONE"="NSAID","IMATINIB"="Kinase inhibitor",
 "VANILLIN"="Other","TADALAFIL"="Other","AZILSARTAN MEDOXOMIL"="Other","LIOTHYRONINE"="Other","SENNOSIDE A"="Other","ESCIN"="Saponin","CYCLOHEXIMIDE"="Protein synthesis inhibitor","MONOBENZONE"="Other")
S$fine <- cls[S$name]
super <- c("Statin"="Sterol / steroid","Corticosteroid"="Sterol / steroid","Sex/steroid hormone"="Sterol / steroid",
 "Antimetabolite"="Antiproliferative","DNA-damaging"="Antiproliferative","Microtubule"="Antiproliferative",
 "Antimicrobial"="Antimicrobial / antiparasitic","Antiparasitic"="Antimicrobial / antiparasitic","CNS/neuroactive"="Neuroactive",
 "Oxidative/thiol"="Oxidative / anti-inflammatory","NSAID"="Oxidative / anti-inflammatory","Kinase inhibitor"="Oxidative / anti-inflammatory",
 "Saponin"="Oxidative / anti-inflammatory","Protein synthesis inhibitor"="Oxidative / anti-inflammatory","Other"="Oxidative / anti-inflammatory")
S$super <- super[S$fine]

## ALL sex-biased hits (no direction filter) -> Venn compartment
S$cat <- ifelse(S$sigWT & !S$sigKD, "WT-only (lost after KD)",
         ifelse(S$sigWT &  S$sigKD, "Shared (Egr1-independent)",
         ifelse(!S$sigWT & S$sigKD, "KD-only (gained after KD)", NA)))
H <- S[!is.na(S$cat) & !is.na(S$super),]
H$cat <- factor(H$cat, levels=c("WT-only (lost after KD)","Shared (Egr1-independent)","KD-only (gained after KD)"))

agg <- as.data.frame(table(H$super, H$cat)); names(agg) <- c("super","cat","n")
tot <- tapply(agg$n, agg$super, sum)
wto <- tapply(agg$n[agg$cat=="WT-only (lost after KD)"], agg$super[agg$cat=="WT-only (lost after KD)"], sum)
pct <- round(100*wto/tot)
ord <- rownames(data.frame(pct,tot))[order(pct, tot)]     # largest % (then larger n) at top
agg$super <- factor(agg$super, levels=ord)
lab <- data.frame(super=factor(names(pct),levels=ord), total=as.integer(tot[names(pct)]),
                  txt=paste0(pct[names(pct)],"% lost after KD"))

th <- theme_bw(base_size=16)+theme(axis.text=element_text(size=14,colour="black"),
      axis.title=element_text(size=16), plot.title=element_text(size=18,face="bold"),
      legend.text=element_text(size=13), legend.title=element_blank(), legend.position="top",
      plot.caption=element_text(size=11,hjust=0,colour="grey30"),
      panel.grid.major.y=element_blank(), panel.grid.minor=element_blank())

p <- ggplot(agg,aes(n,super,fill=cat))+
  geom_col(width=.66,colour="grey30",position=position_stack(reverse=TRUE))+
  geom_text(data=lab,aes(x=total,y=super,label=txt),hjust=-0.08,size=5,fontface="bold",colour="#08519c",inherit.aes=FALSE)+
  scale_fill_manual(values=c("WT-only (lost after KD)"="#1874CD",
                             "Shared (Egr1-independent)"="#a6cee3","KD-only (gained after KD)"="#C71585"))+
  scale_x_continuous(limits=c(0,max(lab$total)+6),expand=c(0,0))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  labs(title="All sex-biased drug hits by class and Egr1 dependence",
       x="Sex-biased drug hits (q ≤ 0.05)",y=NULL,
       caption="Compartments match the panel-C Venn (WT-only 27 | shared 22 | KD-only 3).\nAll hits are male-biased except the 2 KD-only female-biased hits (Guanabenz, Prasterone).")+th
ggsave(file.path(out,"Figure5_panelD_classDependence_allSexBiased.pdf"),p,width=10,height=5.8,device=cairo_pdf)
ggsave(file.path(out,"Figure5_panelD_classDependence_allSexBiased.png"),p,width=10,height=5.8,dpi=200,bg="white")

wide <- reshape(agg,idvar="super",timevar="cat",direction="wide"); wide$total<-tot[as.character(wide$super)]
wide$pct_WTonly <- pct[as.character(wide$super)]
cat("=== all sex-biased hits per class ===\n"); print(wide[order(-wide$pct_WTonly),],row.names=FALSE)
cat("\ncolumn totals -> WT-only:",sum(H$cat=="WT-only (lost after KD)"),
    " Shared:",sum(H$cat=="Shared (Egr1-independent)"),
    " KD-only:",sum(H$cat=="KD-only (gained after KD)")," | grand total:",nrow(H),"\n")
