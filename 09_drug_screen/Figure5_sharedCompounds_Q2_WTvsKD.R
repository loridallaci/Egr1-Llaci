## Figure 5 — the 22 SHARED sex-biased hits: male-vs-female sensitivity (Q2) in
## Egr1 WT vs KD. Shows the sex difference SHRINKS after Egr1 KD while both
## genotypes remain significantly sex-biased (that is what "shared" means).
## Same dumbbell style + Q2 logic as panel f in Figure5_panels_defg.R.
suppressMessages({library(readxl); library(ggplot2); library(dplyr)})
wb  <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Cell Titer Glo/Maxene Drug Screen/85compoundsCheck_analyses_July2026.xlsx"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/09_drug_screen"

d <- suppressMessages(read_excel(wb, sheet="Zc combined", col_names=FALSE))[-1,]
S <- data.frame(name=trimws(as.character(d[[9]])),
                M_WT=as.numeric(d[[3]]), M_KD=as.numeric(d[[4]]),
                F_WT=as.numeric(d[[5]]), F_KD=as.numeric(d[[6]]))
q2 <- function(a,b)(a-b)/sqrt(2)
S$Q2_WT <- q2(S$M_WT,S$F_WT); S$Q2_KD <- q2(S$M_KD,S$F_KD)
S$qWT <- p.adjust(2*pnorm(-abs(S$Q2_WT)),"BH"); S$qKD <- p.adjust(2*pnorm(-abs(S$Q2_KD)),"BH")
S$sigWT <- S$qWT<=0.05; S$sigKD <- S$qKD<=0.05

## SHARED = significant sex-bias in BOTH genotypes (all are male-biased, Q2<0)
sh <- S[S$sigWT & S$sigKD, ]
cat("shared n =", nrow(sh), " | all male-biased:", all(sh$Q2_WT<0), "\n")

## tidy display names (strip salt/ester words, title-case; tag the thioguanine run)
clean <- function(x){
  x <- sub("\\s+(CITRATE|ACETATE|SODIUM|SULFATE|HYDROCHLORIDE|HCL|NITRATE|PHOSPHATE|MESYLATE|PAMOATE|HEMISULFATE|DIPIVOXYL|CALCIUM|PIVALATE|MALEATE|SUCCINATE|POTASSIUM)\\b.*$","",x)
  x <- tools::toTitleCase(tolower(x)); x <- sub("Fluorouracil","Fluorouracil (5-FU)",x); x
}
sh$drug <- clean(sh$name)
sh$drug[sh$name=="THIOGUANINE"] <- "Thioguanine (run 1)"
sh$shift <- sh$Q2_KD - sh$Q2_WT                   # >0 = male bias shrinks toward 0 after KD

## drug-class map (same classes as panel e in Figure5_panels_defg.R)
cls <- c(
 "ATORVASTATIN CALCIUM"="Statin","MEVASTATIN"="Statin","LOVASTATIN"="Statin","SIMVASTATIN"="Statin",
 "PREDNISOLONE SODIUM PHOSPHATE"="Corticosteroid","FLUOCINOLONE ACETONIDE"="Corticosteroid",
 "METHYLPREDNISOLONE"="Corticosteroid","FLUMETHAZONE PIVALATE"="Corticosteroid",
 "FLUDROCORTISONE ACETATE"="Corticosteroid","DEXAMETHASONE"="Corticosteroid",
 "BAZEDOXIFENE ACETATE"="Sex/steroid hormone","TOREMIFENE CITRATE"="Sex/steroid hormone",
 "TAMOXIFEN CITRATE"="Sex/steroid hormone","CLOMIPHENE CITRATE"="Sex/steroid hormone",
 "ESTRONE"="Sex/steroid hormone","ESTRIOL"="Sex/steroid hormone","ESTRADIOL CYPIONATE"="Sex/steroid hormone",
 "ALLYLESTRENOL"="Sex/steroid hormone","NORETHINDRONE ACETATE"="Sex/steroid hormone",
 "PRASTERONE ACETATE"="Sex/steroid hormone","EPLERENONE"="Sex/steroid hormone",
 "ADEFOVIR DIPIVOXYL"="Antimetabolite","AZASERINE"="Antimetabolite","CARMOFUR"="Antimetabolite",
 "FLOXURIDINE"="Antimetabolite","AZACITIDINE"="Antimetabolite","MYCOPHENOLIC ACID"="Antimetabolite",
 "THIOGUANINE"="Antimetabolite","FLUOROURACIL"="Antimetabolite","AZATHIOPRINE"="Antimetabolite",
 "AMSACRINE"="DNA-damaging","OXALIPLATIN"="DNA-damaging","MITOMYCIN"="DNA-damaging","AMINACRINE"="DNA-damaging",
 "DOCETAXEL"="Microtubule","VINBLASTINE SULFATE"="Microtubule",
 "FLUBENDAZOLE"="Antiparasitic","MONENSIN SODIUM"="Antiparasitic","NITARSONE"="Antiparasitic",
 "ARTESUNATE"="Antiparasitic","MEFLOQUINE HYDROCHLORIDE"="Antiparasitic","MEBENDAZOLE"="Antiparasitic",
 "FENBENDAZOLE"="Antiparasitic","PYRVINIUM PAMOATE"="Antiparasitic","ATOVAQUONE"="Antiparasitic",
 "ITRACONAZOLE HYDROCHLORIDE"="Antimicrobial","BENZOXIQUINE"="Antimicrobial","OXICONAZOLE NITRATE"="Antimicrobial",
 "SECNIDAZOLE"="Antimicrobial","DIRITHROMYCIN"="Antimicrobial","OXYQUINOLINE HEMISULFATE"="Antimicrobial",
 "FUSIDIC ACID"="Antimicrobial","FURAZOLIDONE"="Antimicrobial","BENZETHONIUM CHLORIDE"="Antimicrobial",
 "METHYLPHENIDATE HYDROCHLORIDE"="CNS/neuroactive","RILUZOLE"="CNS/neuroactive","METITEPINE MESYLATE"="CNS/neuroactive",
 "KETANSERIN"="CNS/neuroactive","PAROXETINE HYDROCHLORIDE"="CNS/neuroactive","SULOCTIDIL"="CNS/neuroactive",
 "PIMOZIDE"="CNS/neuroactive","FLUPHENAZINE HYDROCHLORIDE"="CNS/neuroactive","GUANABENZ ACETATE"="CNS/neuroactive",
 "ETHOPROPAZINE HYDROCHLORIDE"="CNS/neuroactive",
 "BRONOPOL"="Oxidative/thiol","HYDROQUINONE"="Oxidative/thiol","MENADIONE"="Oxidative/thiol",
 "PHENYLMERCURIC ACETATE"="Oxidative/thiol","THIMEROSAL"="Oxidative/thiol","OLTIPRAZ"="Oxidative/thiol",
 "ASCORBIC ACID"="Oxidative/thiol","CIANIDANOL [+-CATECHIN]"="Oxidative/thiol","ROTENONE"="Oxidative/thiol",
 "CELECOXIB"="NSAID","OXYPHENBUTAZONE"="NSAID","IMATINIB"="Kinase inhibitor",
 "VANILLIN"="Other","TADALAFIL"="Other","AZILSARTAN MEDOXOMIL"="Other","LIOTHYRONINE"="Other",
 "SENNOSIDE A"="Other","ESCIN"="Saponin","CYCLOHEXIMIDE"="Protein synthesis inhibitor","MONOBENZONE"="Other")
sh$class <- cls[sh$name]
if(any(is.na(sh$class))) cat("UNCLASSIFIED:", paste(unique(sh$name[is.na(sh$class)]),collapse=", "),"\n")

## facets ordered top->bottom by how male-biased the class is (mean Q2_WT); drugs within class by Q2_WT
cls_order <- names(sort(tapply(sh$Q2_WT, sh$class, mean)))
sh$class <- factor(sh$class, levels=cls_order)
sh <- sh[order(as.integer(sh$class), -sh$Q2_WT), ]
sh$drug <- factor(sh$drug, levels=sh$drug)

## reduction summary
red <- mean(abs(sh$Q2_WT)) ; redk <- mean(abs(sh$Q2_KD))
nshrunk <- sum(abs(sh$Q2_KD) < abs(sh$Q2_WT))
cat(sprintf("mean |Q2|: WT %.2f -> KD %.2f ; %d/%d compounds shrink toward 0\n",
            red, redk, nshrunk, nrow(sh)))

## glyphs built from code points -> source stays pure ASCII (no encoding issues)
arrL <- intToUtf8(0x2190); arrR <- intToUtf8(0x2192); leq <- intToUtf8(0x2264)

th <- theme_bw(base_size=16)+theme(
  axis.text=element_text(size=14,colour="black"), axis.text.y=element_text(size=14),
  axis.title=element_text(size=16), plot.title=element_text(size=18,face="bold"),
  plot.subtitle=element_text(size=13,colour="grey30",lineheight=1.1),
  axis.title.x=element_text(margin=margin(t=8)),
  legend.text=element_text(size=14), legend.title=element_blank(),
  legend.position="top", panel.grid.minor=element_blank(),
  plot.title.position="plot", plot.margin=margin(10,16,10,10))

p <- ggplot(sh)+
  geom_vline(xintercept=0, colour="grey55")+
  geom_segment(aes(x=Q2_WT, xend=Q2_KD, y=drug, yend=drug),
               arrow=arrow(length=unit(.11,"in")), linewidth=.9, colour="grey55")+
  geom_point(aes(Q2_WT, drug, colour="Egr1 WT"), size=3.6)+
  geom_point(aes(Q2_KD, drug, colour="Egr1 KD"), size=3.6)+
  scale_colour_manual(values=c("Egr1 WT"="#1874CD","Egr1 KD"="grey55"),
                      breaks=c("Egr1 WT","Egr1 KD"),
                      guide=guide_legend(override.aes=list(size=4)))+
  labs(title="Shared hits: sex difference shrinks after Egr1 KD\nbut persists in both genotypes",
       subtitle=sprintf("%d hits significant in both genotypes (q%s0.05)\nmean |Q2| shrinks %.1f %s %.1f (%d/%d compounds) after Egr1 KD",
                        nrow(sh), leq, red, arrR, redk, nshrunk, nrow(sh)),
       x=paste0("Male-vs-female sensitivity (Q2)\n", arrL, " more male-biased          less ", arrR),
       y=NULL)+th+
  facet_grid(rows=vars(class), scales="free_y", space="free_y")+
  theme(strip.text.y=element_text(angle=0, size=12, face="bold", colour="grey15"),
        strip.background=element_rect(fill="grey92", colour="grey70"),
        panel.spacing.y=unit(3,"pt"))

ggsave(file.path(out,"Figure5_sharedCompounds_Q2_WTvsKD.png"), p, width=9, height=11.5, dpi=300, bg="white")
pdf_ok <- tryCatch({ggsave(file.path(out,"Figure5_sharedCompounds_Q2_WTvsKD.pdf"), p, width=9, height=11.5, device=cairo_pdf); TRUE},
                   error=function(e){message("PDF not written (is it open in a viewer?): ", conditionMessage(e)); FALSE})
cat("wrote Figure5_sharedCompounds_Q2_WTvsKD.{pdf,png}\n")
