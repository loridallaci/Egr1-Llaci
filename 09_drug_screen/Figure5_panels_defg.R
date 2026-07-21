## Figure 5 panels d-g, built from verified data:
##  d = Egr1-dependence of the drug sex-bias (49 -> 25 hits after KD)
##  e = drug-class enrichment among male-biased WT hits (hypergeometric vs 85-cpd library)
##  f = statin spotlight: male-vs-female sensitivity (Q2) collapses WT -> KD
##  g = mevalonate/SREBP2 program: log2FC after Egr1 KD, male vs female (+ direct CC targets)
suppressMessages({library(readxl); library(ggplot2); library(dplyr); library(tidyr)})
wb  <- "C:/Users/loril/Documents/Egr1/Egr1 manuscript/Final Submission/Cell Titer Glo/Maxene Drug Screen/85compoundsCheck_analyses_July2026.xlsx"
out <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/09_drug_screen"
de_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/02_rna_analysis/bulk_rna/output_DE_Egr1KD_gRNA3_vs_Neg1"
cc_dir <- "C:/Users/loril/Documents/GitHub/Egr1-Llaci/04_callingcard_analysis/output"
th <- theme_bw(base_size=15)+theme(axis.text=element_text(size=14,colour="black"),
      axis.title=element_text(size=16), plot.title=element_text(size=18,face="bold"),
      legend.text=element_text(size=13), legend.title=element_text(size=14),
      panel.grid.minor=element_blank())

## ---- screen: per-run Q2 in WT and KD ----
d <- suppressMessages(read_excel(wb, sheet="Zc combined", col_names=FALSE))[-1,]
S <- data.frame(name=toupper(trimws(as.character(d[[9]]))),
                M_WT=as.numeric(d[[3]]), M_KD=as.numeric(d[[4]]),
                F_WT=as.numeric(d[[5]]), F_KD=as.numeric(d[[6]]))
q2<-function(a,b)(a-b)/sqrt(2)
S$Q2_WT<-q2(S$M_WT,S$F_WT); S$Q2_KD<-q2(S$M_KD,S$F_KD)
S$qWT<-p.adjust(2*pnorm(-abs(S$Q2_WT)),"BH"); S$qKD<-p.adjust(2*pnorm(-abs(S$Q2_KD)),"BH")
S$sigWT<-S$qWT<=0.05; S$sigKD<-S$qKD<=0.05

## ---- drug -> class map (all 85 runs) ----
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
 "SENNOSIDE A"="Other","ESCIN"="Other","CYCLOHEXIMIDE"="Other","MONOBENZONE"="Other")
S$class <- cls[S$name]
if(any(is.na(S$class))) cat("UNCLASSIFIED:", paste(unique(S$name[is.na(S$class)]),collapse=", "), "\n")

## ================= PANEL d : 49 -> 25 collapse =================
dd <- data.frame(
  State=factor(rep(c("Egr1 WT","Egr1 KD"),each=2),levels=c("Egr1 WT","Egr1 KD")),
  Direction=factor(rep(c("Male-biased","Female-biased"),2),levels=c("Female-biased","Male-biased")),
  n=c(sum(S$sigWT&S$Q2_WT<0),sum(S$sigWT&S$Q2_WT>0),sum(S$sigKD&S$Q2_KD<0),sum(S$sigKD&S$Q2_KD>0)))
pd <- ggplot(dd,aes(State,n,fill=Direction))+geom_col(width=.62,colour="grey25")+
  geom_text(aes(label=ifelse(n>0,n,"")),position=position_stack(vjust=.5),size=6,colour="white",fontface="bold")+
  scale_fill_manual(values=c("Female-biased"="#C71585","Male-biased"="#1874CD"))+
  labs(title="Sex-biased drug response\nis Egr1-dependent",x=NULL,y="Sex-biased hits (q ≤ 0.05)")+th
ggsave(file.path(out,"Figure5_panelD_egr1dependence.pdf"),pd,width=5,height=5.2)
ggsave(file.path(out,"Figure5_panelD_egr1dependence.png"),pd,width=5,height=5.2,dpi=200,bg="white")

## ================= PANEL e : class enrichment among male-biased WT hits =================
hit <- S$sigWT & S$Q2_WT<0
K<-sum(hit); N<-nrow(S)
ce <- S %>% group_by(class) %>% summarise(n=n(), hits=sum(hit[class==class[1]]), .groups="drop")
ce <- S %>% mutate(hit=hit) %>% group_by(class) %>% summarise(n=n(), hits=sum(hit), .groups="drop") %>%
  mutate(frac=hits/n, p=phyper(hits-1,K,N-K,n,lower.tail=FALSE)) %>%
  filter(n>=2) %>% arrange(frac,hits)
ce$class<-factor(ce$class,levels=ce$class)
ce$lab<-sprintf("%d/%d%s",ce$hits,ce$n,ifelse(ce$p<0.05," *",""))
pe <- ggplot(ce,aes(frac,class,fill=-log10(p)))+geom_col(colour="grey30",width=.72)+
  geom_text(aes(label=lab),hjust=-0.12,size=4.6)+
  scale_fill_gradient(low="#deebf7",high="#08519c",name=expression(-log[10]~p))+
  scale_x_continuous(limits=c(0,1.12),breaks=c(0,.5,1),labels=c("0","50%","100%"),expand=c(0,0))+
  labs(title="Drug classes enriched among\nmale-biased hits (WT)",
       x="Fraction of class that is a male-biased hit",y=NULL)+th
ggsave(file.path(out,"Figure5_panelE_classEnrichment.pdf"),pe,width=7.2,height=5.4)
ggsave(file.path(out,"Figure5_panelE_classEnrichment.png"),pe,width=7.2,height=5.4,dpi=200,bg="white")

## ================= PANEL f : statin spotlight =================
st <- c("SIMVASTATIN","LOVASTATIN","MEVASTATIN","ATORVASTATIN CALCIUM")
sf <- S[S$name %in% st,] %>%
  mutate(drug=recode(name,"ATORVASTATIN CALCIUM"="Atorvastatin","SIMVASTATIN"="Simvastatin",
                     "LOVASTATIN"="Lovastatin","MEVASTATIN"="Mevastatin")) %>%
  arrange(Q2_WT)
sf$drug<-factor(sf$drug,levels=sf$drug)
pf <- ggplot(sf)+
  geom_vline(xintercept=0,colour="grey60")+
  geom_segment(aes(x=Q2_WT,xend=Q2_KD,y=drug,yend=drug),arrow=arrow(length=unit(.13,"in")),linewidth=1,colour="grey45")+
  geom_point(aes(Q2_WT,drug,colour="Egr1 WT"),size=4.5)+
  geom_point(aes(Q2_KD,drug,colour="Egr1 KD"),size=4.5)+
  scale_colour_manual(values=c("Egr1 WT"="#1874CD","Egr1 KD"="grey55"),name=NULL,breaks=c("Egr1 WT","Egr1 KD"))+
  labs(title="Statins lose male-selectivity\nafter Egr1 KD",
       x="Male-vs-female sensitivity (Q2)\n← more male-biased",y=NULL)+th
ggsave(file.path(out,"Figure5_panelF_statins.pdf"),pf,width=6,height=4.4)
ggsave(file.path(out,"Figure5_panelF_statins.png"),pf,width=6,height=4.4,dpi=200,bg="white")

## ================= PANEL g : mevalonate/SREBP2 heatmap =================
M<-read.delim(file.path(de_dir,"Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt"))
Fm<-read.delim(file.path(de_dir,"Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"))
ccM<-toupper(na.omit(read.csv(file.path(cc_dir,"overlap_MaleCC_20kb_NearestGene_vs_Egr1g3vsNeg1.csv"))$Gene))
genes<-c("Srebf2","Insig1","Hmgcr","Hmgcs1","Mvk","Mvd","Idi1","Fdps","Fdft1",
         "Sqle","Lss","Cyp51","Dhcr24","Dhcr7","Msmo1","Nsdhl","Sc5d","Ldlr","Vldlr")
gv<-function(sym,tab){r<-tab[!is.na(tab$SYMBOL)&toupper(tab$SYMBOL)==toupper(sym),]
  if(nrow(r)==0)return(c(NA,NA)); r<-r[which.max(r$baseMean),]; c(r$log2FoldChange,r$pvalue)}
H<-do.call(rbind,lapply(genes,function(g){
  vM<-gv(g,M);vF<-gv(g,Fm)
  rbind(data.frame(gene=g,sex="Male",lfc=vM[1],p=vM[2],direct=toupper(g)%in%ccM),
        data.frame(gene=g,sex="Female",lfc=vF[1],p=vF[2],direct=toupper(g)%in%ccM))}))
H$gene<-factor(H$gene,levels=rev(genes)); H$sex<-factor(H$sex,levels=c("Male","Female"))
H$sig<-!is.na(H$p)&H$p<=0.05&abs(H$lfc)>=0.5
ylabs<-rev(sapply(genes,function(g) if(toupper(g)%in%ccM) paste0(g," ★") else g))
pg<-ggplot(H,aes(sex,gene,fill=lfc))+geom_tile(colour="grey85")+
  geom_text(aes(label=ifelse(sig,"*","")),vjust=.78,size=6)+
  scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",midpoint=0,name="log2FC\n(KD vs WT)",limits=c(-1.7,1.7))+
  scale_y_discrete(labels=ylabs)+
  labs(title="Egr1 represses the sterol\nprogram (male-specific)",x=NULL,y=NULL,
       caption="★ direct Egr1 CC target   * DEG (|lfc|≥0.5, p≤0.05)")+
  th+theme(axis.text.y=element_text(size=13),plot.caption=element_text(size=11,hjust=0))
ggsave(file.path(out,"Figure5_panelG_sterolHeatmap.pdf"),pg,width=5,height=7.2)
ggsave(file.path(out,"Figure5_panelG_sterolHeatmap.png"),pg,width=5,height=7.2,dpi=200,bg="white")

cat("built d,e,f,g. WT male-biased hits K =",K,"of N =",N,"\n")
cat("class table:\n"); print(as.data.frame(ce[,c("class","n","hits","frac","p")]),digits=3,row.names=FALSE)
cat("\nstatins:\n"); print(sf[,c("drug","Q2_WT","Q2_KD")],row.names=FALSE,digits=3)
