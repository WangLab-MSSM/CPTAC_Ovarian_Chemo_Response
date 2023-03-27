rm(list=ls())
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(doParallel);nCores=detectCores();registerDoParallel(nCores-2)
library(patchwork)

# setwd("~path")
#
#  ----------- DNA Data Generation  --------------------------------------  
# IDs
meta <- readxl::read_xlsx("Data/Discovery_V2/PTRC_Master_Tumor List_20210212.xlsx", 
                            sheet = "FFPE_Discovery_Samples") %>% 
  janitor::clean_names() 

patient_dna_id <- meta %>%
  select(dna_id, patient_id) %>%
  separate(dna_id, c("dna_id_1", "dna_id_2")) %>% 
  pivot_longer(!patient_id, values_to="dna_id", names_to=NULL) %>%
  filter(dna_id != "NA", !is.na(dna_id))   

patient_rna_id <- meta %>%  select(rna_id, patient_id, dna_id)

# TMB 
mutation <- read_tsv("mutation_all_042922.tsv")

TMB <- mutation %>% 
  gather(patient_id, mutation, -Gene_ID) %>% 
  group_by(patient_id) %>% 
  summarise(tmb = sum(mutation)) %>% 
  mutate(tmb_high = ifelse(tmb >100, 1, 0)) 
TMB %>% count(tmb_high)

# TP53 
TP53= read_csv("tp53_relaxed_mutations.csv") %>% 
  rename(dna_id=case, TP53=tp53_mut)

TP53_patient_level <- TP53 %>%  # take samples that have no duplication
  inner_join(patient_dna_id) %>% 
  group_by(patient_id) %>% 
  summarise(TP53 = 1*any(TP53==T)) 
dim(TP53_patient_level)
# [1] 118   2


# Genome Instability 
GI_table <- readxl::read_xlsx("ptrc-cptac3-oc-cinInfo-it1.xlsx")
tAi5 <- GI_table %>% 
  select(dna_id = sample, tAi5)

tAi5_patient_level <- tAi5 %>%  # take samples that have no duplication
  inner_join(patient_dna_id) %>% 
  group_by(patient_id) %>% 
  summarise(tAi5 = mean(tAi5)) %>%
  select(tAi5, patient_id) %>%
  mutate(tAi5_binary=if_else(tAi5>15, 1, 0))


# chr17LOH 
chr17=read_csv("chr17LOHdata.csv")%>% 
  rename(dna_id="case")

chr17_patient_level <- chr17 %>% 
  inner_join(patient_dna_id) %>% 
  group_by(patient_id) %>% 
  summarise(chr17LOH = max(chr17LOH)) %>% 
  select(patient_id, chr17LOH)

setdiff(TP53_patient_level$patient_id,
        chr17_patient_level$patient_id)
setdiff(chr17_patient_level$patient_id, TP53_patient_level$patient_id)

# BRCA1/2 
brca=read_csv("brca_relaxed_mutations.csv")%>% 
  rename(dna_id="case")
brca_patient_level <- brca %>% 
  inner_join(patient_dna_id) %>% 
  group_by(patient_id) %>% 
  summarise(brca = 1*any(brca_mut==T)) %>% 
  select(patient_id, brca)

# tumor response 
covariates <-  readxl::read_xlsx("PTRC_Master_Tumor List_20210212.xlsx", 
                                 sheet = "FFPE_Discovery_Samples") %>% 
  janitor::clean_names() 

tumor_response <- covariates %>%
  select(patient_id, tumor_response) %>%
  mutate(tumor_response = str_to_sentence(tumor_response))

# put together Genome Instability vs chr17LOH vs TMB 
DNA= TMB  %>% 
  inner_join(tAi5_patient_level) %>% 
  inner_join(chr17_patient_level) %>% 
  inner_join(TP53_patient_level) %>%  
  inner_join(brca_patient_level) %>%
  inner_join(tumor_response) 

write.csv(DNA, file="DNA_patient_level.csv", row.names = F)

#  ----------- Alluvial Plot for DNA events  ------------------------------------------------------

#DNA=read_csv("result/DNA_patient_level.csv")
library(alluvial)

DNA$tmb_100=ifelse(DNA$tmb>100, "High", "Low")
DNA$chr17LOH_b=ifelse(DNA$chr17LOH==1, "Yes", "No")
DNA$TP53_b=ifelse(DNA$TP53==T, "MUT", "1WT")

x=DNA %>% group_by(tumor_response, chr17LOH_b,  brca, TP53_b) %>%
  summarize(n=n())
x2=x[with(x, order(tumor_response, chr17LOH_b, brca,TP53_b)),]
library(ggalluvial)
is_alluvia_form(as.data.frame(x), axes = 1:3, silent = TRUE)


ord <- list(NULL, 
            with(x2, order(tumor_response)), 
            with(x2, order(tumor_response)), 
            with(x2, order(tumor_response)))
pdf("Paper_PTRC/plot/Fig2F_right.pdf",width = 4,height = 3)
alluvial(select(x2,   chr17LOH_b, tumor_response, brca, TP53_b), 
         col = ifelse(x$tumor_response == "Sensitive",  "#FFC20A", "#0C7BDC"),
         freq=x2$n, ordering=ord, cw=0.15)
dev.off();



# Regression Data generation  -----------------------------------------------------------------
# cov: age, S/R, tumor location, purity
PTRC_Purity_estimates <- read_tsv("PTRC_Purity_estimates.tsv") %>% 
  separate(sample, c("a", "b", "c", "d", "e"),"_") %>% 
  select(d,ESTIMATE) %>% 
  mutate(d = str_replace_all(d, "\\.", "-")) %>% 
  rename(patient_id = d) %>% 
  group_by(patient_id) %>% 
  summarise(purity = mean(ESTIMATE))

covariates <- readxl::read_xlsx("/PTRC_Master_Tumor List_20201201.xlsx", sheet = "FFPE_Discovery_Samples") %>% 
  janitor::clean_names() %>% 
  select(patient_id = patient_id, patient_age, tumor_location_group, tumor_response) %>% 
  mutate(tumor_location_group = case_when(tumor_location_group == "OM" ~ "OM",
                                          tumor_location_group == "OV" ~ "OV",
                                          TRUE ~ "Others")) %>% 
  mutate(tumor_response = str_to_lower(tumor_response)) %>% 
  inner_join(PTRC_Purity_estimates)
# protein & rna 
     # Take mean for duplicated samples
protein_regression <- read_tsv("protein_regression_iprofun.tsv")
protein_regression_long <- protein_regression %>% 
  gather(patient_id, protein, -Gene_ID)

rna_regression <- read_tsv("rna_regression_iprofun.tsv")
rna_regression_long <- rna_regression %>% 
  gather(patient_id, rna, -Gene_ID)

dat_protein_DNA <- protein_regression_long %>% 
  inner_join(DNA %>% select(!tumor_response)) %>% 
  inner_join(covariates) 
GeneID=unique(dat_protein_DNA$Gene_ID)

dat_protein_rna_DNA <- rna_regression_long  %>% 
  full_join(protein_regression_long)  %>% 
  inner_join(DNA%>% select(!tumor_response)) %>% 
  left_join(covariates) 

save(DNA, dat_protein_DNA, dat_protein_rna_DNA, GeneID, file="DNA_cov_rna_protein.RData")



#  ----------- Protein/Pathway Association Analysis (not stratified by S/R) ------------------------------------------------------
# regression DNA (CHR17LOH, BRCA) vs proteins  
#load("Paper_PTRC/Rdata/DNA_cov_rna_protein.RData")
a <- dat_protein_DNA %>% 
  group_by(Gene_ID) %>% 
  nest() %>% 
  mutate(model = map(data,~lm(protein ~ tmb_high + chr17LOH + brca+ patient_age + 
                                purity + tumor_location_group+ tumor_response, data = .)))  %>% 
  mutate(model_info = map(model, broom::tidy)) %>%
  dplyr::select(-data, -model) %>%
  unnest %>% 
  ungroup


reg.summary=a %>% filter(term=="chr17LOH") %>%
  mutate(FDR_BH = p.adjust(p.value, method = "BH")) %>%
  write_csv("chr17LOH_adj_protein_3cat.csv")

reg.summary=a %>% filter(term=="brca") %>%
  mutate(FDR_BH = p.adjust(p.value, method = "BH")) %>%
  write_csv("brca_adj_protein_3cat.csv")


# pathway analysis with regression summaries 
source("code/function_pathway.R") 

chr17LOH.summary=read_csv("Paper_PTRC/result/chr17LOH_adj_protein_3cat.csv")
chr17LOH.pathway=logp_pathway_multiple(GeneID=chr17LOH.summary$Gene_ID, 
                        EST=chr17LOH.summary$estimate, 
                        Pvalue=chr17LOH.summary$p.value)  %>%
  write_csv("chr17LOH.adj.pathway.csv")

brca.summary=read_csv("Paper_PTRC/result/brca_adj_protein_3cat.csv")
brca.pathway=logp_pathway_multiple(GeneID=brca.summary$Gene_ID, 
                                       EST=brca.summary$estimate, 
                                       Pvalue=brca.summary$p.value)  %>%
  write_csv("brca.adj.pathway.csv")


pathway_protein <- chr17LOH.pathway  %>%select(geneset.names, chr17LOH_avglogFC=avglogFC,  
                                        chr17LOH_fdr=fdr.BY) %>%
  full_join(brca.pathway  %>%select(geneset.names, BRCA_avglogFC=avglogFC,  
                                   BRCA_fdr=fdr.BY)) %>%
  write_csv("pathway_tran_protein.csv")


#  ----------- TP53WT score - Protein/Pathway Association Analysis (not stratified by S/R) ------------------------------------------------------
# regression TP53WT score vs proteins  
load("Rdata/DNA_cov_rna_protein.RData")
TP53score=read.csv("Data/TP53_activity_scores.csv")

dat2=TP53score %>% full_join(patient_rna_id, by=c("RNA_ID"="rna_id")) %>%
  group_by(patient_id) %>% summarise(TP53score.wt=mean(ssGSEA.TP53.WT, na.rm = T))%>% 
  mutate(TP53WT.score.high=if_else(TP53score.wt>0.25, 1, 0)) %>%
  select(patient_id, TP53WT.score.high) %>% inner_join(dat_protein_DNA)


a <- dat2 %>% 
  group_by(Gene_ID) %>% 
  nest() %>% 
  mutate(model = map(data,~lm(protein ~ TP53WT.score.high+ patient_age + 
                                purity + tumor_location_group+ tumor_response, data = .)))  %>% 
  mutate(model_info = map(model, broom::tidy)) %>%
  dplyr::select(-data, -model) %>%
  unnest %>% 
  ungroup

reg.summary=a %>% filter(term=="TP53WT.score.high") %>%
  mutate(FDR_BH = p.adjust(p.value, method = "BH")) %>%
  write_csv("Revision/Result/TP53WT.score.high_adj_protein.csv")


# pathway analysis with regression summaries 
source("code/function_pathway.R") 

tp53.summary=read_csv("Revision/Result/TP53WT.score.high_adj_protein.csv")
tp53.pathway=logp_pathway_multiple(GeneID=tp53.summary$Gene_ID, 
                                  EST=tp53.summary$estimate, 
                                  Pvalue=tp53.summary$p.value) %>%
  write_csv("Revision/Result/tp53.adj.pathway.csv")




# stratified analysis    ----------------------------------------------------------------- 
load("Rdata/DNA_cov_rna_protein.RData")

idx=c("chr17LOH", "brca")
TP53score=read.csv("Data/TP53_activity_scores.csv")

dat2=TP53score %>% inner_join(patient_rna_id, by=c("RNA_ID"="rna_id")) %>%
  group_by(patient_id) %>% summarise(TP53score.wt=mean(ssGSEA.TP53.WT, na.rm = T))%>% 
  mutate(TP53WT.score.high=if_else(TP53score.wt>0.25, 1, 0)) %>%
  select(patient_id, TP53WT.score.high) %>% full_join(dat_protein_DNA)



# S/R vs protein stratified by TMB
StratifiedTrans=function(strata) {
  Coef=foreach (i = 1:length(GeneID), .combine="rbind") %dopar% {
    dat=dat2 %>% filter(Gene_ID == GeneID[i])
    
    dat$Low_RvsS=(dat[,strata]==0)*(dat$tumor_response=="refractory")
    dat$High_RvsS=(dat[,strata]==1)*(dat$tumor_response=="refractory")
    temp=summary(lm(protein~ Low_RvsS + 
                      High_RvsS  + 
                      patient_age + purity +  tumor_location_group, data=dat))$coefficients[2:3,]
    return(c(temp[1,], temp[2,]))
  }
  colnames(Coef)=c("Low.est", "Low.se", "Low.T", "Low.pvalue", "High.est", "High.se", "High.T", "High.pvalue")
  Coef=data.frame(Coef)
  Coef2=cbind(GeneID, Coef)
  Coef2$Low.BH=p.adjust(Coef2$Low.pvalue, "BH")
  Coef2$High.BH=p.adjust(Coef2$High.pvalue, "BH")
  Coef3=Coef2[, c("GeneID", "Low.est", "Low.se", "Low.T", "Low.pvalue", "Low.BH", 
                  "High.est", "High.se", "High.T", "High.pvalue", "High.BH")]
  write.csv(Coef3, 
            file=paste0("Revision/result/Protein_vs_SR_stratified_", strata, ".csv"), 
            row.names=F)
  
}

StratifiedTrans("chr17LOH")
StratifiedTrans("brca")
StratifiedTrans("TP53WT.score.high")



# Stratified plot of L1CAM  -----------------------------------------------------------------
load("Paper_PTRC/Rdata/DNA_cov_rna_protein.RData")
res <- read_csv(paste0("Protein_vs_SR_stratified_chr17LOH.csv"))%>%
  filter(GeneID=="L1CAM") %>% select(Low.pvalue, High.pvalue) %>% 
  as.matrix() %>% format(., scientific=T,digits=2)
dat=dat_protein_rna_DNA %>% 
  filter(Gene_ID == "L1CAM") %>% 
  mutate(tumor_response=factor(tumor_response, levels=c("sensitive", "refractory"),
                               labels=c("Sensitive", "Refractory"))) 
p_protein_L1CAM <- dat %>% 
  ggplot(aes(x = as.factor(chr17LOH), y = protein, color =tumor_response ))+
  geom_boxplot()+
  labs(x = "chr17-LOH", y = "L1CAM Protein Level")+
  theme_pubr()+ 
  theme(legend.position = "none",)+
  scale_color_manual(name="", values=c("#FFC20A", "#0C7BDC"))+ 
  annotate(geom="text", x=1, y=3.4, label= res[1])+ 
  annotate(geom="text", x=2, y=3.4, label=res[2])+ 
  scale_x_discrete(labels=c("0" = "No", "1" = "Yes"))+
  ggbeeswarm::geom_quasirandom(dodge.width=.8) 
p_protein_L1CAM

# 
pdf(file="Paper_PTRC/plot/L1CAM.pdf", 
    width =2.25, height = 4, paper="a4r")
p_protein_L1CAM
dev.off()




# Stratified plot of TGM2  -----------------------------------------------------------------

# TGM2
load("Paper_PTRC/Rdata/DNA_cov_rna_protein.RData")

res <- read_csv(paste0("Revision/Result/Protein_vs_SR_stratified_TP53WT.score.high.csv")) %>%
  filter(GeneID=="TGM2") %>% select(Low.pvalue, High.pvalue) %>% 
  as.matrix() %>% format(., scientific=T,digits=2)

TP53score=read.csv("TP53_activity_scores.csv")
dat1=dat_protein_DNA%>% 
  filter(Gene_ID == "TGM2") %>% 
  mutate(tumor_response=factor(tumor_response, levels=c("sensitive", "refractory"),
                               labels=c("Sensitive", "Refractory")))


dat2=TP53score %>% inner_join(patient_rna_id, by=c("RNA_ID"="rna_id")) %>%
  group_by(patient_id) %>% summarise(TP53score.wt=mean(ssGSEA.TP53.WT, 
                                                       na.rm = T))%>% 
  mutate(TP53WT.score.high=if_else(TP53score.wt>0.25, 1, 0)) %>%
  select(patient_id, TP53score.wt, TP53WT.score.high) %>% inner_join(
    dat1)
p_TGM2 <- dat2 %>% 
  ggplot(aes(x = as.factor(TP53WT.score.high), y = protein, color =
               tumor_response))+
  geom_boxplot()+
  theme_pubr()+ 
  scale_x_discrete(labels=c("1" = "High", "0" = "Low"))+
  theme(legend.position = "none") +
  scale_color_manual(values=c("#FFC20A", "#0C7BDC"))+ 
  labs(x = "TP53-WT score", y= "TGM2 relative protein level") + 
  annotate(geom="text", x=1, y=3.4, label= res[1])+ 
  annotate(geom="text", x=2, y=3.4, label=res[2])+
  ggbeeswarm::geom_quasirandom(dodge.width=.8) 
p_TGM2


# TGM2

pdf(file="Revision/Result/TGM2_plot.pdf", 
    width =2.5, height = 4, paper="a4r")
p_TGM2
dev.off()

#  ----------- Plot Pathways Stratified Analysis 


idx=c("chr17LOH", "TP53WT.score.high")


  dat <- read_csv(paste0("Revision/Result/Protein_vs_SR_stratified_", idx[i], ".csv"))  %>%
    mutate(label=case_when(Low.BH <0.1| High.BH<0.1 ~ GeneID))
  dat$Direction="1Not_Sig"
  dat$Direction[which(dat$Low.est>0 & dat$High.est>0 & (dat$Low.BH <0.1| dat$High.BH<0.1))]="High_Positive_Low_Positive"
  dat$Direction[which(dat$Low.est<0 & dat$High.est>0 & (dat$Low.BH <0.1| dat$High.BH<0.1))]="High_Positive_Low_Negative"
  dat$Direction[which(dat$Low.est>0 & dat$High.est<0 & (dat$Low.BH <0.1| dat$High.BH<0.1))]="High_Negative_Low_Positive"
  dat$Direction[which(dat$Low.est<0 & dat$High.est<0 & (dat$Low.BH <0.1| dat$High.BH<0.1))]="High_Negative_Low_Negative"
  
  color_pool=cbind(c("gray", "red",  "orange", "purple", "#0072B2"), 
                   c("1Not_Sig", "High_Positive_Low_Positive", 
                     "High_Positive_Low_Negative", "High_Negative_Low_Positive", 
                     "High_Negative_Low_Negative"))
  color_use=color_pool[match(names(table(dat$Direction)), color_pool[,2]),1]
  
  
  # TP53WT.score.high
  pic=
    ggplot(dat, aes(x=-log10(Low.pvalue), 
                    y=-log10(High.pvalue), 
                    label=GeneID, 
                    color=Direction)) + 
    geom_point() + theme_pubr()+
    scale_color_manual(values =  color_use,
                       name = "Association", 
                       labels = c("Not Significant", 
                                  "Neg(Low)/Neg(High)", 
                                  "Pos(Low)/Neg(High)",
                                  "Pos(Low)/Pos(High)"))+
    ggrepel::geom_text_repel(aes(label = label)) +
    # geom_text(aes(label=label,
    #           position=position_jitter(), hjust=0,vjust=0)+
    scale_x_continuous(limits = c(00,6))+
    scale_y_continuous(limits = c(00,6))+
    geom_abline(intercept=0, slope=1, col=2)+ 
    labs(x = "-log10(p) in low TP53-WT score",
         y= "-log10(p) in high TP53-WT score")+
    theme(legend.justification=c(1,1),legend.position=c(1,1),
          axis.title.x = element_text(size = 15), 
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15), 
          axis.text.y = element_text(size = 15))
  pdf(file="Revision/Result/Protein_vs_SR_stratified_TP53WT.score.pdf", 
      width =4, height = 4, paper="a4r")
  pic
  dev.off()
  
  # chr17LOH
  pic=
    ggplot(dat, aes(x=-log10(Low.pvalue), 
                    y=-log10(High.pvalue), 
                    label=GeneID, 
                    color=Direction)) + 
    geom_point() + theme_pubr() +
    scale_color_manual(values =  color_use,
                       name = "Association", 
                       labels = c("Not Significant", 
                                  "Neg/Neg", 
                                  "Pos/Pos"))+
    ggrepel::geom_text_repel(aes(label = label)) +
    scale_x_continuous(limits = c(00,6))+
    scale_y_continuous(limits = c(00,6))+
    geom_abline(intercept=0, slope=1, col=2)+ 
    labs(x = "-log10(p) in no chr17-LOH ",
         y= "-log10(p) in chr17-LOH")+
    theme(axis.title.x = element_text(size = 15), 
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          legend.justification=c(1,1),legend.position=c(1,1))
  pdf(file="Paper_PTRC/plot/Protein_vs_SR_chr17LOH.pdf", 
      width =4, height = 4, paper = "a4r")
  pic  
  dev.off()
  pic
  







