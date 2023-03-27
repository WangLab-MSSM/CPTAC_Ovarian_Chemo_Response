rm(list=ls())
library(data.table)
library(tidyverse)
library(iProFun)

# setwd("~path")
# load data ---------------------------------------------------------------
# narrow analysis to those with group info.
group=read.table("CPTAC2_clustering.txt", header=T)
temp1=matrix(unlist(strsplit(group$Sample, "_")), ncol=2,byrow = T)
temp2=matrix(unlist(strsplit(temp1[,2], "[.]")), ncol=3,byrow = T)
group$SubID=SubID=paste0("TCGA-", temp2[,1], "-", temp2[,2])

# outcome
rna=fread("retrospective_ova_array_sort_common_gene_15121.txt")
rna_subject=rna[,intersect(c("Gene_ID", SubID), colnames(rna)), with=F] %>% as.data.frame()

protein=fread("retro_ova_LRcombined_proteome_sort_common_gene_7061.txt")
protein_subject=protein[,intersect(c("Gene_ID", SubID), colnames(protein)), with=F] %>% as.data.frame()

# cnv 
cnv = fread("retrospective_ova_CNA_sort_common_gene_11859.txt") 
cnv_subject=cnv[,intersect(c("Gene_ID", SubID), colnames(cnv)), with=F] %>% as.data.frame()

cnv_subject[,group[group$Response=="Refractory", SubID]]
cnv_subject_refractory=cnv_subject_sensitive= cnv_subject
cnv_subject_sensitive[,group[group$Response=="Refractory", "SubID"]]=0
cnv_subject_refractory[,group[group$Response=="Sensitive", "SubID"]]=0
# methy
# create gene-level methylation data from probe level
methy=fread("retrospective_OV_methy27_cleaned.txt")
NAs=apply(methy[,c(-1:-3)], 1, function(f) mean(1*is.na(f)))
methy2=methy[which(NAs<0.9),]
methy_prob_subject=methy2[,intersect(c("Gene_ID", SubID), colnames(methy2)), with=F] %>% as.data.frame()
methy_subject1=sapply(2:ncol(methy_prob_subject),
                  function(f) tapply(methy_prob_subject[,f], methy_prob_subject$Gene_ID, mean))
methy_subject=methy_subject1 %>% as.data.frame()%>%rownames_to_column(var="Gene_ID")
colnames(methy_subject)=colnames(methy_prob_subject)
methy_subject_refractory=methy_subject_sensitive= methy_subject
methy_subject_sensitive[,group[group$Response=="Refractory", "SubID"]]=0
methy_subject_refractory[,group[group$Response=="Sensitive", "SubID"]]=0

# cov - NULL
cov_sensitive=group %>% 
  mutate(Sensitive=ifelse(Response=="Sensitive", 1, 0)) %>%
         select("SubID", "Sensitive")   %>%
  column_to_rownames("SubID") %>% t()%>% as.data.frame()


cov_refractory=group %>% 
  mutate(Refractory=ifelse(Response=="Refractory", 1, 0)) %>%
  select("SubID", "Refractory")   %>%
  column_to_rownames("SubID") %>% t()%>% as.data.frame()





#  iProFun input with refractory as reference ---------------------------------------------------------------

yList = list(
  rna_subject  %>% as.data.frame() ,
  protein_subject %>% as.data.frame()
)
xList = list(cnv_subject %>% as.data.frame(),
             methy_subject %>% as.data.frame(),
             cnv_subject_sensitive %>% as.data.frame(),
             methy_subject_sensitive %>% as.data.frame()
)
covariates = list(
  cov_sensitive %>% as.data.frame(),
  cov_sensitive %>% as.data.frame()
)

pi1 =0.05 

# Reg
reg.refractory.ref=iProFun.reg(yList=yList, xList, covariates, 
                               var.ID=c("Gene_ID"), Y.rescale=T)


# Calculate Empirical FDR for all outcomes
eFDR.refractory.ref=iProFun.eFDR(reg.all=reg.refractory.ref, yList=yList, xList=xList, 
                                 covariates=covariates, pi1=pi1,
                                 NoProbXIndex=NULL, Y.rescale = T,
                                 permutate_number=5, var.ID=c("Gene_ID"), seed=123)
# iProFun identification

res.refractory.ref=iProFun.detection(reg.all=reg.refractory.ref, 
                                     eFDR.all=eFDR.refractory.ref, 
                                     FWER.all=NULL, 
                                     filter=c(1, -1, 0, 0),
                                     NoProbButFWERIndex=NULL, 
                                     fdr.cutoff = 0.1, fwer.cutoff=0.1, 
                                     PostPob.cutoff=0.75,
                                     xType=c("cnv_ref", "methy_ref", "cnv_inter", 
                                             "methy_inter"), 
                                     yType=c("rna", "protein"))


#  iProFun input with sensitive as reference ---------------------------------------------------------------

yList = list(
  rna_subject  %>% as.data.frame() ,
  protein_subject %>% as.data.frame()
)
xList = list(cnv_subject %>% as.data.frame(),
             methy_subject %>% as.data.frame(),
             cnv_subject_refractory %>% as.data.frame(),
             methy_subject_refractory %>% as.data.frame()
)
covariates = list(
  cov_refractory %>% as.data.frame(),
  cov_refractory %>% as.data.frame()
)

pi1 =0.05 

# Reg
reg.sensitive.ref=iProFun.reg(yList=yList, xList, covariates, 
                               var.ID=c("Gene_ID"), Y.rescale=T)


# Calculate Empirical FDR for all outcomes
eFDR.sensitive.ref=iProFun.eFDR(reg.all=reg.sensitive.ref, yList=yList, xList=xList, 
                                 covariates=covariates, pi1=pi1,
                                 NoProbXIndex=NULL, Y.rescale = T,
                                 permutate_number=5, var.ID=c("Gene_ID"), seed=126)

reg.sum=multi.omic.reg.summary(reg.out.list=reg.sensitive.ref)
reg.tab=iProFun.reg.table(reg.all=reg.sensitive.ref, 
                          xType = c("cnv_ref", "methy_ref", "cnv_inter", "methy_inter"), 
                          yType = c("rna", "protein"))

# iProFun identification

res.sensitive.ref=iProFun.detection(reg.all=reg.sensitive.ref, 
                                     eFDR.all=eFDR.sensitive.ref, 
                                     FWER.all=NULL, 
                                     filter=c(1, -1, 0, 0),
                                     NoProbButFWERIndex=NULL, 
                                     fdr.cutoff = 0.1, fwer.cutoff=0.1, 
                                     PostPob.cutoff=0.75,
                                     xType=c("cnv_ref", "methy_ref", "cnv_inter", 
                                             "methy_inter"), 
                                     yType=c("rna", "protein"))


# harvest results  ------------------------------------------------------------------

result = res.refractory.ref %>% 
  mutate(Tissue=if_else(xType=="cnv_ref" | xType=="methy_ref", 
                        "Refractory", "Interaction")) %>%
  filter(Tissue=="Refractory") %>% 
  bind_rows(res.sensitive.ref %>% 
              mutate(Tissue=if_else(xType=="cnv_ref" | xType=="methy_ref", 
                                    "Sensitive", "Interaction"))
  ) %>%
  mutate(xType=case_when(xType=="cnv_ref" ~ "cnv",
                         xType=="methy_ref" ~ "methy",
                         xType=="cnv_inter" ~ "cnv",
                         xType=="methy_inter" ~ "methy"))

write_tsv(result, file="CPTAC2_iProFun_result_021023.tsv")
save(result, file= "CPTAC2_iProFun_result_021023.Rdata")



# iPrFun pathway plots   ------------------------------------------------------
library(tidyverse)
library(iProFun)
library(epitools)
load("Paper_PTRC/Revision/Result/CPTAC2_iProFun_result_021023.Rdata")
dat1=result  %>% filter(xType=="cnv") %>% filter(yType=="rna" | yType=="protein")
idx0=dat1 %>%  filter(Tissue =="Refractory") %>% select(xName) 
Genes_2platform=data.frame(idx0[duplicated(idx0),] )[,1]
#  4795 (PTRC 8584) genes with both mRNA and protein 

idx1 = dat1  %>% filter(iProFun.identification==1  & Tissue =="Refractory") %>% 
  select(xName) 
Genes_R_cascade = data.frame(idx1[duplicated(idx1),] )[,1]
Genes_R_not = setdiff(Genes_2platform, Genes_R_cascade)
idx2 = dat1  %>% filter(iProFun.identification==1  & Tissue =="Sensitive") %>% 
  select(xName) 
Genes_S_cascade = data.frame(idx2[duplicated(idx2),])[,1] 
Genes_S_not = setdiff(Genes_2platform, Genes_S_cascade)
# Sensitive: 2348 (PTRC 4682)
# Refractory : 1948 (PTRC 3606)
pr=length(Genes_R_cascade)/(length(Genes_R_not)+length(Genes_R_cascade))
# 0.4062565 (PTRC 0.4200839) 
ps=length(Genes_S_cascade)/(length(Genes_S_not)+length(Genes_S_cascade))
# 0.4896767 (PTRC 0.5454334)
res <- prop.test(x = c(2348, 1948), n = c(4795, 4795))

# Pathways 
source("Paper_PTRC/code/function_pathway.R")
R_result=binary_pathway_multiple(YesList=Genes_R_cascade, NoList=Genes_R_not) 
S_result=binary_pathway_multiple(YesList=Genes_S_cascade, NoList=Genes_S_not) 
iProFun_pathway=cbind(S_result[,1:3], S=S_result[,4:8], R=R_result[,4:8])  
iProFun_pathway$geneset.names=gsub(" ", "_", iProFun_pathway$geneset.names )


iProFun_pathway2 = data.frame(geneset.names="All",
                              n_path=length(Genes_2platform),
                              n_path_data=length(Genes_2platform),
                              S.prop.yes.in.path=ps,
                              R.prop.yes.in.path=pr
                              ) %>% 
  bind_rows(iProFun_pathway)
 

write_csv(iProFun_pathway2, "CPTAC2_cascade_genes_pathway_021023.csv")


# ------------- validation pathway plot-----------------------------

library(ggpubr)
pic=path2 %>% right_join(format) %>%
  ggplot(aes(x = s_prop_yes_in_path, y = r_prop_yes_in_path)) +
  geom_point(aes(size =n_path_data, color=as.factor(bubble_color), alpha = alpha))+
  theme_pubr()+
  ggrepel::geom_text_repel(aes(label = pathway_label , color =as.factor(text_color)),
                           show.legend  = FALSE, size = 3, position=position_jitter()) +
  geom_abline( slope = 1, intercept = 0, color = "grey") +
  
  scale_size(range = c(1, 20), breaks=c(100, 300))+ 
  scale_color_manual(values =  c("gray25","red", "blue", "black" , "#E68613"), guide = 'none')+
  
  
  scale_x_continuous(limits = c(0,0.8))+
  scale_y_continuous(limits = c(0,0.7))+
  guides(fill=FALSE)+
  labs(x = "% CNV-RNA/Protein cascades in CPTAC2 proxy-sensitive tumor",
       y = "% CNV-RNA/Protein cascades in CPTAC2 proxy-refractory tumor",
       size= "Pathway Size")  + theme(legend.position = c(0.9, 0.115))+ 
  scale_alpha(guide = 'none')
pic

pdf("CPTAC2_iProFun_pathway.pdf",width = 6.3,height =5,paper="a4r")
pic  
dev.off()









