rm(list=ls())

library(tidyverse)
library(iProFun)
library(epitools)

# setwd("~path")
#---------------------------------------------------------------#
# ---------------- Data cleaning for iProFun Analysis ------------------------------
#---------------------------------------------------------------#
# ------------------ id matching ------------------#

meta <- readxl::read_xlsx("PTRC_Master_Tumor List_20210212.xlsx", 
                              sheet = "FFPE_Discovery_Samples") %>% janitor::clean_names() 

patient_dna_id <- meta %>%  select(dna_id, patient_id) %>%
  separate(dna_id, c("dna_id_1", "dna_id_2")) %>% 
  pivot_longer(!patient_id, values_to="dna_id", names_to=NULL) %>%
  filter(dna_id != "NA", !is.na(dna_id))   


patient_rna_id <- meta %>%  select(rna_id, patient_id) %>%
  separate(rna_id, c("rna_id_1", "rna_id_2")) %>% 
  pivot_longer(!patient_id, values_to="rna_id", names_to=NULL) %>%
  filter(rna_id != "NA", !is.na(rna_id))  

meta.pro <- readxl::read_xlsx("PTRC_Master_Tumor List_20210212.xlsx", 
                          sheet = "FD_TMT_order") %>% janitor::clean_names() 
patient_prot_id <- meta.pro %>% select(patient_id, 
                                      pro_id=data_label_2_data_set_tm_tplex_sample_source_patient_id_sample_id)
  
# ---------- mutation --------------------------------- #

mutation_strict <- read_csv("mutation_all-proteins-strict.csv") %>% 
  janitor::clean_names() %>% 
  select(Gene_ID = symbol,
         dna_id) %>% 
  mutate(value = 1) %>% 
  drop_na()

mutation_new <- readxl::read_xlsx("CPTAC3-OC-recurrently-mutated-hgsc-pathogenic-and-cosmic-hotspots.xlsx") %>% 
  janitor::clean_names() %>% 
  select(Gene_ID = symbol,
         dna_id) %>% 
  mutate(value = 1) %>% 
  drop_na() 


mutation_combined <- mutation_strict %>% 
  filter(!(Gene_ID %in% intersect(unique(mutation_new$Gene_ID),
                                unique(mutation_strict$Gene_ID)))) %>% 
  bind_rows(mutation_new)

mutation_all <- mutation_combined %>% 
  unique %>% 
  left_join(patient_dna_id) %>% 
  drop_na() %>% 
  group_by(Gene_ID, patient_id) %>% 
  summarise(value = max(value)) %>% 
  ungroup() %>% 
  spread(patient_id, value, fill = 0) %>% 
  write_tsv("mutation_all_042922.tsv")

mutation=mutation_all %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  filter(sum >= 6) %>% 
  select(-sum) %>% 
  write_tsv("mutation_iprofun_042922.tsv")

# rbind two raw tables for data release

mutation_new <- readxl::read_xlsx("CPTAC3-OC-recurrently-mutated-hgsc-pathogenic-and-cosmic-hotspots.xlsx") %>% 
  janitor::clean_names() %>% 
  mutate(source="AddManualStep") %>% 
  drop_na(dna_id) %>% 
  mutate(aft=as.numeric(aft), adt=as.numeric(adt), dpt=as.numeric(dpt), adt_fwd=as.numeric(adt_fwd),
         adt_rev=as.numeric(adt_rev), str=as.numeric(str))

GeneManual=unique(mutation_new$symbol)

mutation_strict <- read_csv("mutation_all-proteins-strict.csv") %>% 
  janitor::clean_names() %>%
  mutate(source="Computation") %>%
  filter(!symbol %in% GeneManual) %>%
  select(!c("afn", "adn", "dpn"))

mutation_combined=mutation_new  %>% bind_rows(mutation_strict) %>% 
  write_tsv("mutation_fulltable_071322.tsv")

# ---------- cnv -------------------------------------- #
cnv <- read_tsv("FD_CNV_Log2Freq_20210104_v02.tsv")  %>% 
  select(Gene_ID = `HGNC symbol`, starts_with("D")) %>% 
  pivot_longer(!Gene_ID, names_to="dna_id", values_to = "cnv") %>%
  unique %>% left_join(patient_dna_id) %>% 
  drop_na() %>% group_by(Gene_ID, patient_id) %>% 
  summarise(cnv = mean(cnv)) %>% 
  ungroup() %>% 
  spread(patient_id, cnv) %>% 
  write_tsv("cnv_iprofun_042922.tsv")
sum(is.na(cnv))
dim(cnv) # 19920   121
cnv <- cnv[apply(cnv, 1, function(f) mean(is.na(f))<=0.5),] 
dim(cnv) # 19920   121

# --------------------------------- loh ------------------------------ #
loh <- read_csv("loh_v4.1.csv") %>% 
  rename(Gene_ID=gene_name) %>%
  pivot_longer(!Gene_ID, names_to="dna_id") %>%
  unique %>% left_join(patient_dna_id) %>% 
  drop_na() %>% group_by(Gene_ID, patient_id) %>% 
  summarise(value = max(value)) %>% 
  ungroup() %>% 
  spread(patient_id, value) %>% 
  write_tsv("loh_iprofun_042922.tsv")
sum(is.na(loh))
dim(loh) # 32598   121
loh <- loh[apply(loh, 1, function(f) mean(is.na(f))<=0.5),] 
dim(loh) # 32598   121
 
# ------------------------------ rna --------------------------------- # 
rna <- read_tsv("FD_RNAseq_Log2cpm_20210215_v02.tsv")  %>% 
  select(Gene_ID = `HGNC symbol`, starts_with("R")) %>%
  pivot_longer(!Gene_ID, names_to="rna_id") %>%
  unique %>% left_join(patient_rna_id) %>% 
  drop_na() %>% group_by(Gene_ID, patient_id) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  spread(patient_id, value) %>% 
  write_tsv("rna_iprofun_042922.tsv")

dim(rna) # 19605   107
rna <- rna[apply(rna, 1, function(f) mean(is.na(f))<=0.5),] 
dim(rna)# 19605   107

# ----------------------------- protein --------------------------------- # 
protein <- read_tsv("FD_GLBLprot_MI_FDbridge_Abund_20201002_Imput_v02.tsv") %>% 
  select(Gene_ID = `HGNC symbol`, starts_with("FD")) %>%
  pivot_longer(!Gene_ID, names_to="pro_id") %>%
  unique %>% left_join(patient_prot_id) %>% 
  drop_na() %>% group_by(Gene_ID, patient_id) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  spread(patient_id, value) %>% 
  write_tsv("protein_iprofun_042922.tsv")

# --------------------------- phospho ----------------------------------- # 
phospho <- read_tsv("FD_PHOSsite_MI_FDbridge_MsiteAbund_20201030_Imput_v02.tsv") %>% 
  select(Gene_ID = `HGNC symbol`, Site=`Index`, starts_with("FD")) %>%
  pivot_longer(cols=-c(Gene_ID, Site), names_to="pro_id") %>%
  unique %>% left_join(patient_prot_id) %>% 
  drop_na() %>% group_by(Gene_ID, Site, patient_id) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  spread(patient_id, value) %>% 
  write_tsv("phospho_iprofun_042922.tsv")
  
# ---------------------------  covariates -----------------------------------  #  
PTRC_Purity_estimates <- read_tsv("PTRC_Purity_estimates.tsv") %>% 
  separate(sample, c("a", "b", "c", "d", "e"),"_") %>% 
  select(d,ESTIMATE) %>% 
  mutate(d = str_replace_all(d, "\\.", "-")) %>% 
  rename(patient_id = d) %>% 
  group_by(patient_id) %>% 
  summarise(purity = mean(ESTIMATE))


qc_warning_rna <- readxl::read_xlsx("cptac3-rna-quality-control.xlsx") %>% 
  janitor::clean_names() %>% 
  rename(rna_id = library) %>% 
  inner_join(patient_rna_id) %>% 
  select(patient_id, qc_warning) %>% 
  mutate(qc_warning = case_when(is.na(qc_warning) ~ 0,
                                TRUE ~ 1)) 

cov.pro.pool <- meta %>%
  select(patient_id, patient_age, tumor_response, tumor_type) %>%
  mutate(tumor_response = case_when(tumor_response == "Refractory" ~ "refractory",
                                    tumor_response == "Sensitive" ~ "sensitive",
                                    TRUE ~ tumor_response)) %>%
  mutate(tumor_type = case_when(tumor_type == "Metastatic" ~ 2,
                                tumor_type == "Primary" ~ 0,
                                tumor_type == "Mix" ~ 1,
  )) %>% 
  mutate(refractory = ifelse(tumor_response == "refractory", 1, 0)) %>% 
  mutate(sensitive = ifelse(tumor_response == "sensitive", 1, 0))%>% 
  select(-tumor_response) %>% 
  inner_join(PTRC_Purity_estimates) 

cov.rna.pool = cov.pro.pool%>%  inner_join(qc_warning_rna)

cov.pro.pool.wide = cov.pro.pool%>%
  pivot_longer(!patient_id, names_to = "key", values_to = "value") %>%
  spread(key = "patient_id", value = "value") %>% 
  column_to_rownames(var = "key") 
cov.rna.pool.wide = cov.rna.pool%>%
  pivot_longer(!patient_id, names_to = "key", values_to = "value") %>%
  spread(key = "patient_id", value = "value") %>% 
  column_to_rownames(var = "key") 

cov.pro.sensitive.ref=  cov.pro.pool.wide%>% filter(
  row.names(cov.pro.pool.wide) %in% c("patient_age","tumor_type", 
                                  "refractory", "purity"))  %>% 
  write_tsv("covariates_iprofun_pro_sensitive_ref_042922.tsv")

cov.rna.sensitive.ref  = cov.rna.pool.wide  %>% filter(
    row.names(cov.rna.pool.wide) %in% c("patient_age","tumor_type", 
                               "refractory", "purity", "qc_warning"))  %>%   
  write_tsv("covariates_iprofun_rna_sensitive_ref_042922.tsv")

cov.rna.refractory.ref= cov.rna.pool.wide %>% filter(
  row.names(cov.rna.pool.wide) %in% c("patient_age","tumor_type", 
                             "sensitive", "purity", "qc_warning"))  %>%   
  write_tsv("covariates_iprofun_rna_refractory_ref_042922.tsv")

cov.pro.refractory.ref= cov.pro.pool.wide %>% filter(
  row.names(cov.pro.pool.wide) %in% c("patient_age","tumor_type", 
                             "sensitive", "purity"))  %>% 
  write_tsv("covariates_iprofun_pro_refractory_ref_042922.tsv")


# ------------------ DNA interaction  --------------------------------- # 
sensitive = cov.pro.pool %>% select(patient_id, sensitive)
refractory = cov.pro.pool %>% select(patient_id, refractory)

loh.sensitive <- loh %>% 
  pivot_longer(!Gene_ID, names_to="patient_id") %>%
  inner_join(sensitive) %>% 
  mutate(interaction=value * sensitive) %>%
  select(Gene_ID, patient_id, interaction) %>%
  spread(key = "patient_id", value = "interaction") %>% 
  write_tsv("loh_sensitive_iprofun_042922.tsv")

cnv.sensitive <- cnv %>% 
  pivot_longer(!Gene_ID, names_to="patient_id") %>%
  inner_join(sensitive) %>% 
  mutate(interaction=value * sensitive) %>%
  select(Gene_ID, patient_id, interaction) %>%
  spread(key = "patient_id", value = "interaction") %>% 
  write_tsv("cnv_sensitive_iprofun_042922.tsv")

mutation.sensitive <- mutation %>% 
  pivot_longer(!Gene_ID, names_to="patient_id") %>%
  inner_join(sensitive) %>% 
  mutate(interaction=value * sensitive) %>%
  select(Gene_ID, patient_id, interaction) %>%
  spread(key = "patient_id", value = "interaction") %>% 
  write_tsv("mutation_sensitive_iprofun_042922.tsv")


loh.refractory <- loh %>% 
  pivot_longer(!Gene_ID, names_to="patient_id") %>%
  inner_join(refractory) %>% 
  mutate(interaction=value * refractory) %>%
  select(Gene_ID, patient_id, interaction) %>%
  spread(key = "patient_id", value = "interaction") %>% 
  write_tsv("loh_refractory_iprofun_042922.tsv")

cnv.refractory <- cnv %>% 
  pivot_longer(!Gene_ID, names_to="patient_id") %>%
  inner_join(refractory) %>% 
  mutate(interaction=value * refractory) %>%
  select(Gene_ID, patient_id, interaction) %>%
  spread(key = "patient_id", value = "interaction") %>% 
  write_tsv("cnv_refractory_iprofun_042922.tsv")

mutation.refractory <- mutation %>% 
  pivot_longer(!Gene_ID, names_to="patient_id") %>%
  inner_join(refractory) %>% 
  mutate(interaction=value * refractory) %>%
  select(Gene_ID, patient_id, interaction) %>%
  spread(key = "patient_id", value = "interaction") %>% 
  write_tsv("mutation_refractory_iprofun_042922.tsv")


save(mutation, loh, cnv, 
     mutation.sensitive, cnv.sensitive, loh.sensitive,
     mutation.refractory, cnv.refractory, loh.refractory,
     cov.pro.pool,
     cov.rna.sensitive.ref, cov.pro.sensitive.ref, 
     cov.rna.refractory.ref, cov.pro.refractory.ref,
     rna, protein, phospho,
     file="iProFun_data.rda")

#---------------------------------------------------------------#
# ---------------- iProFun Analysis ------------------------------
#---------------------------------------------------------------#

#  ------------ iProFun input with refractory as reference --------------------#

yList = list(
  rna  %>% as.data.frame() ,
  protein %>% as.data.frame(),
  phospho  %>% as.data.frame()
)
xList = list(cnv %>% as.data.frame(),
             loh %>% as.data.frame(),
             mutation %>% as.data.frame(),
             cnv.sensitive %>% as.data.frame(),
             loh.sensitive %>% as.data.frame(),
             mutation.sensitive %>% as.data.frame()
)
covariates = list(
  cov.rna.refractory.ref %>% as.data.frame(),
  cov.pro.refractory.ref %>% as.data.frame(),
  cov.pro.refractory.ref %>% as.data.frame()
)

pi1 =0.05 

# Reg
reg.refractory.ref=iProFun.reg(yList=yList, xList, covariates, 
                               var.ID=c("Gene_ID"), var.ID.additional = c("Site"), Y.rescale=T)


#  FWER for data type(s) that have few number of genes
FWER.refractory.ref=iProFun.FWER(reg.all=reg.refractory.ref, FWER.Index=c(3, 6))

# Calculate Empirical FDR for all outcomes
eFDR.refractory.ref=iProFun.eFDR(reg.all=reg.refractory.ref, yList=yList, xList=xList, 
                                 covariates=covariates, pi1=pi1,
                                 NoProbXIndex=c(3,6), Y.rescale = T,
                                 permutate_number=1, var.ID=c("Gene_ID"),
                                 var.ID.additional=c( "Site"), seed=123)
# iProFun identification

res.refractory.ref=iProFun.detection(reg.all=reg.refractory.ref, eFDR.all=eFDR.refractory.ref, 
                                     FWER.all=FWER.refractory.ref, filter=c(1, 0, 0, 0, 0, 0),
                                     NoProbButFWERIndex=c(3,6), fdr.cutoff = 0.1, fwer.cutoff=0.1, 
                                     PostPob.cutoff=0.75,
                                     xType=c("cnv", "loh", "mutation", "cnv_sensitive", 
                                             "loh_sensitive", "mutation_sensitive"), 
                                     yType=c("rna", "protein", "phospho"))

# repeat the analysis with a different reference (sensitive = 0)
xList = list(cnv %>% as.data.frame(),
             loh %>% as.data.frame(),
             mutation %>% as.data.frame(),
             cnv.refractory %>% as.data.frame(),
             loh.refractory %>% as.data.frame(),
             mutation.refractory %>% as.data.frame()
)
covariates = list(
  cov.rna.sensitive.ref %>% as.data.frame(),
  cov.pro.sensitive.ref %>% as.data.frame(),
  cov.pro.sensitive.ref %>% as.data.frame()
)


# Reg
reg.sensitive.ref=iProFun.reg(yList=yList, xList, covariates, 
                              var.ID=c("Gene_ID"), var.ID.additional = c("Site"), Y.rescale=T)


#  FWER for data type(s) that have few number of genes
FWER.sensitive.ref=iProFun.FWER(reg.all=reg.sensitive.ref, FWER.Index=c(3, 6))

# Calculate Empirical FDR for all outcomes
eFDR.sensitive.ref=iProFun.eFDR(reg.all=reg.sensitive.ref, yList=yList, xList=xList, 
                                covariates=covariates, pi1=pi1,
                                NoProbXIndex=c(3,6), Y.rescale = T,
                                permutate_number=1, var.ID=c("Gene_ID"),
                                var.ID.additional=c( "Site"), seed=123)
# iProFun identification

res.sensitive.ref=iProFun.detection(reg.all=reg.sensitive.ref, eFDR.all=eFDR.sensitive.ref, 
                                    FWER.all=FWER.sensitive.ref, filter=c(1, 0, 0, 0, 0, 0),
                                    NoProbButFWERIndex=c(3,6), fdr.cutoff = 0.1, fwer.cutoff=0.1, 
                                    PostPob.cutoff=0.75,
                                    xType=c("cnv", "loh", "mutation", "cnv_sensitive", 
                                            "loh_sensitive", "mutation_sensitive"), 
                                    yType=c("rna", "protein", "phospho"))


# --------------------- harvest results  --------------------------- # 

result = res.refractory.ref %>% 
  mutate(Tissue=if_else(xType=="mutation" | xType=="cnv" | xType=="loh", 
                        "Refractory", "Interaction")) %>%
  filter(Tissue=="Refractory") %>% 
  bind_rows(res.sensitive.ref %>% 
              mutate(Tissue=if_else(xType=="mutation" | xType=="cnv" | xType=="loh", 
                                    "Sensitive", "Interaction"))
  ) %>%
  mutate(xType=case_when(xType=="cnv_sensitive" ~ "cnv",
                         xType=="loh_sensitive" ~ "loh",
                         xType=="mutation_sensitive" ~ "mutation",
                         TRUE ~ xType))

write_tsv(result, file="Paper_PTRC/result/latest iProFun result/iProFun_result_042922.tsv")
save(result, file= "iProFun_result_042922.Rdata")
