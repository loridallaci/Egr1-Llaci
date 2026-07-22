## CSV version of the "Sex-biased drug-sensitivity hits by class and Egr1 dependence" table.
## Compartments = q<=0.05 Male-vs-Female sensitivity hits (Q2=(M-F)/sqrt(2), BH-FDR over 85 runs):
##   WT-only = sig in Egr1 WT only (lost after KD) ; Shared = sig in both ; KD-only = sig in KD only (gained).
suppressMessages(library(readxl))
wb  <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Cell Titer Glo/Maxene Drug Screen/85compoundsCheck_analyses_July2026.xlsx"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/09_drug_screen"

d <- suppressMessages(read_excel(wb, sheet="Zc combined", col_names=FALSE))[-1,]
S <- data.frame(name=toupper(trimws(as.character(d[[9]]))),
                M_WT=as.numeric(d[[3]]), M_KD=as.numeric(d[[4]]),
                F_WT=as.numeric(d[[5]]), F_KD=as.numeric(d[[6]]))
q2 <- function(a,b)(a-b)/sqrt(2)
S$Q2_WT <- q2(S$M_WT,S$F_WT); S$Q2_KD <- q2(S$M_KD,S$F_KD)
S$q_WT  <- p.adjust(2*pnorm(-abs(S$Q2_WT)),"BH"); S$q_KD <- p.adjust(2*pnorm(-abs(S$Q2_KD)),"BH")
S$sigWT <- S$q_WT<=0.05; S$sigKD <- S$q_KD<=0.05

## drug-class map (same fine classes as the analysis)
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

## compartment + direction
S$Compartment <- ifelse(S$sigWT & !S$sigKD, "WT-only (lost after KD)",
                 ifelse(S$sigWT &  S$sigKD, "Shared (Egr1-independent)",
                 ifelse(!S$sigWT & S$sigKD, "KD-only (gained after KD)", NA)))
hits <- S[!is.na(S$Compartment), ]
## direction from the significant genotype (WT if WT-sig, else KD)
hits$Direction <- ifelse(ifelse(hits$sigWT, hits$Q2_WT, hits$Q2_KD) < 0, "Male-biased", "Female-biased")
hits$Class <- cls[hits$name]

## tidy display name; tag the two thioguanine runs in order of appearance
clean <- function(x){
  x <- sub("\\s+(CITRATE|ACETATE|SODIUM|SULFATE|HYDROCHLORIDE|HCL|NITRATE|PHOSPHATE|MESYLATE|PAMOATE|HEMISULFATE|DIPIVOXYL|CALCIUM|PIVALATE|MALEATE|SUCCINATE|POTASSIUM|MEDOXOMIL)\\b.*$","",x)
  tools::toTitleCase(tolower(x)) }
hits$Drug <- clean(hits$name)
hits$Drug[hits$name=="FLUOROURACIL"] <- "Fluorouracil (5-FU)"
tio <- which(hits$name=="THIOGUANINE"); if(length(tio)==2) hits$Drug[tio] <- paste0("Thioguanine (run ", seq_along(tio), ")")

hits <- hits[order(factor(hits$Compartment,
          levels=c("WT-only (lost after KD)","Shared (Egr1-independent)","KD-only (gained after KD)")),
          hits$Class, hits$Q2_WT), ]
csv <- hits[, c("Drug","Class","Compartment","Direction","Q2_WT","q_WT","Q2_KD","q_KD")]
csv[,c("Q2_WT","q_WT","Q2_KD","q_KD")] <- lapply(csv[,c("Q2_WT","q_WT","Q2_KD","q_KD")], function(v) signif(v,3))
outfile <- file.path(out,"Figure5_drugClass_table.csv")
write.csv(csv, outfile, row.names=FALSE)

cat("wrote", outfile, "\n")
cat("compartment counts:\n"); print(table(hits$Compartment))
cat("total hits:", nrow(hits), " | male-biased:", sum(hits$Direction=="Male-biased"),
    " female-biased:", sum(hits$Direction=="Female-biased"), "\n")
