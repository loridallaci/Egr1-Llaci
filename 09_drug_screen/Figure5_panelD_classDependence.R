## Figure 5 panel D — "Sex-biased drug vulnerability is broadly Egr1-dependent,
## except antiproliferatives". Stacked bar of male-biased sex-difference hits per
## drug super-class, split by Egr1-dependence. Bars ordered by % Egr1-dependent (desc).
##
## HOW THE NUMBERS ARE MADE:
##  - Source = 85-run confirmation screen ("Zc combined" sheet), M/F x Egr1-WT/KD averaged Z.
##  - Sex-difference per genotype: Q2 = (Z_male - Z_female)/sqrt(2); two-sided normal p; BH-FDR over 85 runs.
##  - A run is a "male-biased hit" in a genotype if q<=0.05 AND Q2<0 (males more sensitive).
##  - Each hit is assigned to ONE Egr1-dependence category:
##       lost after KD  = male-biased in WT but NOT in KD  (Egr1-dependent)
##       shared         = male-biased in BOTH WT and KD    (Egr1-independent)
##       gained after KD= male-biased in KD but NOT in WT
##  - Drugs are grouped into super-classes; % Egr1-dep = lost / (all male-biased hits in that class).
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

## super-class map; the old catch-all "Other" is renamed to a descriptive label
super <- c("Statin"="Sterol / steroid","Corticosteroid"="Sterol / steroid","Sex/steroid hormone"="Sterol / steroid",
 "Antimetabolite"="Antiproliferative","DNA-damaging"="Antiproliferative","Microtubule"="Antiproliferative",
 "Antimicrobial"="Antimicrobial / antiparasitic","Antiparasitic"="Antimicrobial / antiparasitic",
 "CNS/neuroactive"="Neuroactive",
 "Oxidative/thiol"="Oxidative / anti-inflammatory","NSAID"="Oxidative / anti-inflammatory","Kinase inhibitor"="Oxidative / anti-inflammatory",
 "Saponin"="Oxidative / anti-inflammatory","Protein synthesis inhibitor"="Oxidative / anti-inflammatory","Other"="Oxidative / anti-inflammatory")
S$super <- super[S$fine]

## male-biased hits and their Egr1-dependence category
S$cat <- ifelse(S$sigWT&S$Q2_WT<0 & !(S$sigKD&S$Q2_KD<0), "Egr1-dependent (lost after KD)",
         ifelse(S$sigWT&S$Q2_WT<0 &  (S$sigKD&S$Q2_KD<0), "Egr1-independent (shared)",
         ifelse(!(S$sigWT&S$Q2_WT<0) & S$sigKD&S$Q2_KD<0,  "Gained after KD", NA)))
H <- S[!is.na(S$cat),]
H$cat <- factor(H$cat, levels=c("Egr1-dependent (lost after KD)","Egr1-independent (shared)","Gained after KD"))

## per-super counts + % Egr1-dependent, ordered by % (largest at top)
agg <- as.data.frame(table(H$super, H$cat)); names(agg) <- c("super","cat","n")
tot <- tapply(agg$n, agg$super, sum)
dep <- tapply(agg$n[agg$cat=="Egr1-dependent (lost after KD)"], agg$super[agg$cat=="Egr1-dependent (lost after KD)"], sum)
pct <- round(100*dep/tot)
ord <- names(sort(pct))                              # ascending -> largest ends up at TOP of the plot
agg$super <- factor(agg$super, levels=ord)
lab <- data.frame(super=factor(names(pct),levels=ord), total=as.integer(tot[names(pct)]),
                  txt=paste0(pct[names(pct)],"% Egr1-dep"))

th <- theme_bw(base_size=16)+theme(axis.text=element_text(size=14,colour="black"),
      axis.title=element_text(size=16), plot.title=element_text(size=18,face="bold"),
      legend.text=element_text(size=13), legend.title=element_blank(), legend.position="top",
      panel.grid.major.y=element_blank(), panel.grid.minor=element_blank())

p <- ggplot(agg,aes(n,super,fill=cat))+
  geom_col(width=.66,colour="grey30")+
  geom_text(data=lab,aes(x=total,y=super,label=txt),hjust=-0.08,size=5,fontface="bold",colour="#08519c",inherit.aes=FALSE)+
  scale_fill_manual(values=c("Egr1-dependent (lost after KD)"="#1874CD",
                             "Egr1-independent (shared)"="#a6cee3","Gained after KD"="#C71585"))+
  scale_x_continuous(limits=c(0,max(lab$total)+5),expand=c(0,0))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  labs(title="Sex-biased drug vulnerability is broadly Egr1-dependent,\nexcept antiproliferatives",
       x="Male-biased drug hits",y=NULL)+th
ggsave(file.path(out,"Figure5_panelD_classDependence.pdf"),p,width=10,height=5.5,device=cairo_pdf)
ggsave(file.path(out,"Figure5_panelD_classDependence.png"),p,width=10,height=5.5,dpi=200,bg="white")

cat("=== panel D counts ===\n")
wide <- reshape(agg,idvar="super",timevar="cat",direction="wide"); wide$total<-tot[as.character(wide$super)]
wide$pct <- pct[as.character(wide$super)]; print(wide[order(-wide$pct),], row.names=FALSE)
