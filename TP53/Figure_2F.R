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

#  ----------- Alluvial Plot for DNA events : 2F ------------------------------------------------------

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

