## Drug-class table collapsed to the 5 super-classes used in panel D, x Egr1-dependence
## compartment (WT-only / Shared / KD-only). Same idea as Figure5_drugClass_table.pdf.
suppressMessages({library(readxl); library(gridExtra); library(grid)})
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
super <- c("Statin"="Sterol / steroid","Corticosteroid"="Sterol / steroid","Sex/steroid hormone"="Sterol / steroid",
 "Antimetabolite"="Antiproliferative","DNA-damaging"="Antiproliferative","Microtubule"="Antiproliferative",
 "Antimicrobial"="Antimicrobial / antiparasitic","Antiparasitic"="Antimicrobial / antiparasitic","CNS/neuroactive"="Neuroactive",
 "Oxidative/thiol"="Oxidative / anti-inflammatory","NSAID"="Oxidative / anti-inflammatory","Kinase inhibitor"="Oxidative / anti-inflammatory",
 "Saponin"="Oxidative / anti-inflammatory","Protein synthesis inhibitor"="Oxidative / anti-inflammatory","Other"="Oxidative / anti-inflammatory")
S$super <- super[cls[S$name]]

S$comp <- ifelse(S$sigWT & !S$sigKD, "WT",
          ifelse(S$sigWT &  S$sigKD, "Shared",
          ifelse(!S$sigWT & S$sigKD, "KD", NA)))
H <- S[!is.na(S$comp) & !is.na(S$super),]

## display names
clean <- function(x){ x<-sub("\\s+(CITRATE|ACETATE|SODIUM|SULFATE|HYDROCHLORIDE|NITRATE|PHOSPHATE|MESYLATE|PAMOATE|HEMISULFATE|CALCIUM|PIVALATE|MALEATE|MEDOXOMIL|ACETONIDE)\\b.*$","",x)
  tools::toTitleCase(tolower(x)) }
H$Drug <- clean(H$name)
H$Drug[H$name=="FLUOROURACIL"]        <- "Fluorouracil (5-FU)"
H$Drug[H$name=="PRASTERONE ACETATE"]  <- "Prasterone (DHEA)"
H$Drug[H$name=="PHENYLMERCURIC ACETATE"] <- "Phenylmercuric acetate"
tio <- which(H$name=="THIOGUANINE"); if(length(tio)==2) H$Drug[tio] <- paste0("Thioguanine (run ", seq_along(tio), ")")

super_order <- c("Sterol / steroid","Oxidative / anti-inflammatory","Antimicrobial / antiparasitic","Neuroactive","Antiproliferative")
cell <- function(sc, cc){ x<-H[H$super==sc & H$comp==cc,]
  if(nrow(x)==0) return("—")
  key <- if(cc=="KD") x$Q2_KD else x$Q2_WT
  paste(x$Drug[order(key)], collapse="\n") }
mat <- sapply(c("WT","Shared","KD"), function(cc) sapply(super_order, function(sc) cell(sc,cc)))
df <- data.frame(cls=super_order, mat, check.names=FALSE, stringsAsFactors=FALSE)
colnames(df) <- c("Drug class",
  "WT-only  (n = 27)\nlost after Egr1 KD","Shared  (n = 22)\nEgr1-independent","KD-only  (n = 3)\ngained after KD")

nr <- nrow(df); nc <- ncol(df)
colfill  <- c("#f0f0f0","#EAF3FB","#F4FAFD","#EAF6F9")            # body column tints
headfill <- c("#d9d9d9","#6BAED6","#C6DBEF","#3690C0")            # header fills
tt <- ttheme_minimal(
  core   = list(fg_params=list(hjust=0, x=0.04, vjust=1, y=0.94, fontsize=12, col="grey10"),
                bg_params=list(fill=rep(colfill, each=nr), col="grey85", lwd=1)),
  colhead= list(fg_params=list(fontsize=13, fontface="bold", col=c("black","black","black","white")),
                bg_params=list(fill=headfill, col="grey60", lwd=1)))
g <- tableGrob(df, rows=NULL, theme=tt)
## bold the class-name column
lay <- g$layout
for(r in 2:(nr+1)){ id <- which(lay$t==r & lay$l==1 & lay$name=="core-fg"); if(length(id)) g$grobs[[id]]$gp <- gpar(fontface="bold", fontsize=12) }
## sensible column widths
g$widths <- unit(c(2.6, 3.1, 3.1, 3.1), "in")

title <- textGrob("Sex-biased drug-sensitivity hits by class and Egr1 dependence",
                  gp=gpar(fontsize=18, fontface="bold"))
foot  <- textGrob(paste0("Compartments = q ≤ 0.05 Male-vs-Female sensitivity hits (BH-FDR over 85 runs), grouped into the 5 panel-D classes.\n",
                         "All hits are male-biased except the 2 KD-only female-biased hits (Guanabenz, Prasterone). Thioguanine was screened twice."),
                  gp=gpar(fontsize=11, col="grey30"), just="left", x=0.02, hjust=0)
tabh <- sum(as.numeric(grid::convertHeight(g$heights,"in")))
G <- arrangeGrob(title, g, foot, ncol=1, heights=unit(c(0.5, tabh, 0.75), "in"))
W <- sum(as.numeric(grid::convertWidth(g$widths,"in")))+0.3; Hh <- 0.5+tabh+0.75+0.3

cairo_pdf(file.path(out,"Figure5_drugClass_table_5groups.pdf"), width=W, height=Hh); grid.draw(G); dev.off()
png(file.path(out,"Figure5_drugClass_table_5groups.png"), width=W, height=Hh, units="in", res=200, bg="white"); grid.draw(G); dev.off()

cat("wrote Figure5_drugClass_table_5groups.{pdf,png}\n")
tab <- table(factor(H$super,levels=super_order), factor(H$comp,c("WT","Shared","KD")))
print(addmargins(tab))
