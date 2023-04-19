output<-read.csv("FFPE_Protein_Assoc.tsv",sep="\t") #### change to FFPE rna assoc, Frozen protein assoc, FFPE validation protein/rna assoc, CPTAC-2016 protein/rna assoc

global.vol<-output
global<-global.vol
genep<-cbind.data.frame(global.vol$symbol,global.vol$p.val_sen)

max(genep[,2],na.rm=T)
###################################

library(qusage)

go<-read.gmt(file="pathway_March1/Hallmark.gmt") #### change to KEGG.gmt, Reactome.gmt, DDR.gmt, DQ.gmt, TGFb.gmt

geneID<-genep[,1]
# --- map genes to terms
go.final<-matrix(0,length(geneID),length(go))
colnames(go.final)<-names(go); rownames(go.final)<-geneID

for (j in 1:length(go)) {
  m.g<-which(geneID %in% go[[j]])
  go.final[m.g,j]<-1
}

pathway<-go.final

sum<-NULL
for(i in 1:ncol(pathway))
{
  s=sum(pathway[,i])
  sum<-c(sum,s)
}
ind<-which(sum>4)

sub.pathway<-pathway[,ind]
dim(sub.pathway)

#################################

mean.vec.hitps<-NULL
mean.vec.missps<-NULL
pval.vec<-NULL
size<-NULL

for(i in 1:dim(sub.pathway)[2])
{
  print(i)
  ind<-which(sub.pathway[,i]==1)
  genes<-rownames(sub.pathway)[ind]
  size<-c(size,length(genes))
  
  hitp<--log(global$p.val_sen[which(global$HGNC.symbol %in% genes)])/log(10)
  hitf<-global$FC_sen[which(global$HGNC.symbol %in% genes)]
  signhitp<-ifelse(hitf>1,1,-1)
  hitps<-hitp*signhitp
  mean.vec.hitps<-c(mean.vec.hitps,mean(hitps,na.rm=T))
  
  missp<--log(global$p.val_sen[-(which(global$HGNC.symbol %in% genes))])/log(10)
  missf<-global$FC_sen[-(which(global$HGNC.symbol %in% genes))]
  signmissp<-ifelse(missf>1,1,-1)
  missps<-missp*signmissp
  mean.vec.missps<-c(mean.vec.missps,mean(missps,na.rm=T))
  
  wst<-wilcox.test(hitps,missps,alternative = "two.sided",paired=FALSE,mu=0)
  
  pval.wst<-wst$p.value
  pval.vec<-c(pval.vec,pval.wst)
}

updown<-ifelse(mean.vec.hitps>mean.vec.missps,"up","down")
fdr.wst=p.adjust(pval.vec,method="BH",n=length(pval.vec))

signed_pathway<-cbind.data.frame("pathway"=paste0("TGFb_",colnames(sub.pathway)),"size"=size,"pvalue"=pval.vec,"fdr"=fdr.wst,"mean_in"=mean.vec.hitps,"mean_out"=mean.vec.missps,"updown"=updown)

write.table(signed_pathway,file="signed_pathway_Hallmark.tsv",sep="\t",quote = F,row.names = F, col.names = T)

############ phospho collapsing site to gene ############
phos<-read.csv("FFPE_phospho_assoc.tsv",sep="\t") #### change to Frozen/FFPE validation/CPTAC-2016 phospho

phos$symbol<- sub("\\..*", "", as.character(as.factor(phos$symbol)))

colnames(phos)[2]<-"GeneID"
l=length(unique(phos$GeneID))
u=unique(phos$GeneID)

pval.pos<-NULL
pval.neg<-NULL

coef.pos<-NULL
coef.neg<-NULL

for(i in 1:l)
{
  print(i)
  ind<-which(phos$GeneID %in% u[i])
  sub.phos<-phos[ind,]
  ind.pos<-which(sub.phos$FC_sen>1)
  ind.neg<-which(sub.phos$FC_sen<1)
  
  if(length(ind.pos)>0)
  {
    sub.phos.pos<-sub.phos[ind.pos,]
    
    m<-which.min(sub.phos.pos$p.val_sen)
    
    pval.pos<-c(pval.pos,sub.phos.pos$p.val_sen[m])
    coef.pos<-c(coef.pos,sub.phos.pos$FC_sen[m])
  }else{
    pval.pos<-c(pval.pos,1)
    coef.pos<-c(coef.pos,NA)
  }
  if(length(ind.neg)>0)
  {
    sub.phos.neg<-sub.phos[ind.neg,]
    
    m<-which.min(sub.phos.neg$p.val_sen)
    
    pval.neg<-c(pval.neg,sub.phos.neg$p.val_sen[m])
    coef.neg<-c(coef.neg,sub.phos.neg$FC_sen[m])
  }else{
    pval.neg<-c(pval.neg,1)
    coef.neg<-c(coef.neg,NA)
  }
}

asso.pos<-cbind.data.frame(u,pval.pos,coef.pos)
colnames(asso.pos)<-c("Gene","p.val","FC")

asso.neg<-cbind.data.frame(u,pval.neg,coef.neg)
colnames(asso.neg)<-c("Gene","p.val","FC")

write.table(asso.pos,file="phos_pos.tsv",sep="\t",quote = F,row.names = F, col.names = T)
write.table(asso.neg,file="phos_neg.tsv",sep="\t",quote = F,row.names = F, col.names = T)


######## load "phos_pos.tsv" ########
pval.vec<-NULL
stat.vec<-NULL
size<-NULL

for(i in 1:dim(sub.pathway)[2])
{
  print(i)
  ind<-which(sub.pathway[,i]==1)
  genes<-rownames(sub.pathway)[ind]
  size<-c(size,length(genes))
  
  hitp<-global$p.val[which(global$Gene %in% genes)]
  
  missp<-global$p.val[-(which(global$Gene %in% genes))]

  wst<-wilcox.test(hitp,missp,alternative = "less",paired=FALSE,mu=0)
  
  stat.wst<-wst$statistic
  stat.vec<-c(stat.vec,stat.wst)
  
  pval.wst<-wst$p.value
  pval.vec<-c(pval.vec,pval.wst)
}

fdr.wst=p.adjust(pval.vec,method="BH",n=length(pval.vec))

pathway.pos<-cbind.data.frame("pathway"=colnames(sub.pathway),"size.pos"=size,"pvalue.pos"=pval.vec,"fdr.pos"=fdr.wst)

####### load "phos_neg.tsv" and run the enrichment test again ############
pathway.neg<-cbind.data.frame("pathway"=colnames(sub.pathway),"size.neg"=size,"pvalue.neg"=pval.vec,"fdr.neg"=fdr.wst)

pval.fin<-ifelse(pathway.pos$fdr.pos < pathway.neg$fdr.neg,pathway.pos$pvalue.pos,pathway.neg$pvalue.neg)
fdr.fin<-ifelse(pathway.pos$fdr.pos < pathway.neg$fdr.neg,pathway.pos$fdr.pos,pathway.neg$fdr.neg)  
size.fin<-ifelse(pathway.pos$fdr.pos < pathway.neg$fdr.neg,pathway.pos$size.pos,pathway.neg$size.neg)   
dir.fin<-ifelse(pathway.pos$fdr.pos < pathway.neg$fdr.neg,"up","down")  

pathway.fin<-cbind.data.frame("pathway"=colnames(sub.pathway),"size"=size.fin,"pvalue"=pval.fin,"fdr"=fdr.fin,"updown"=dir.fin)


write.table(pathway.fin,file="signed_pathway_Hallmark.tsv",sep="\t",quote = F,row.names = F, col.names = T)


########## pathway dotplot ############

prot<- read.csv("FFPE_global_pathway.tsv", sep="\t", header=T)
phos<- read.csv("FFPE_phospho_pathway.tsv", sep="\t", header=T)
rna<- read.csv("FFPE_RNA_pathway.tsv", sep="\t", header=T)
path<- read.csv("pathway_marcin.csv", header=T)
path1<- path[!is.na(path$short),]

prot.sub<- prot[prot$pathway %in% path1$pathway,]
prot.sub.m<- merge(prot.sub,path1, by="pathway")
prot.sub.m$short<- paste0(prot.sub.m$short," (",prot.sub.m$Total.genes,")")

phos.sub<- phos[phos$pathway %in% path1$pathway,]
phos.sub.m<- merge(phos.sub,path1, by="pathway")
phos.sub.m$short<- paste0(phos.sub.m$short," (",phos.sub.m$Total.genes,")")

rna.sub<- rna[rna$pathway %in% path1$pathway,]
rna.sub.m<- merge(rna.sub,path1, by="pathway")
rna.sub.m$short<- paste0(rna.sub.m$short," (",rna.sub.m$Total.genes,")")

sig_prot<-cbind.data.frame(prot.sub.m$pathway,-log(prot.sub.m$fdr),prot.sub.m$updown,prot.sub.m$fdr, prot.sub.m$short)
colnames(sig_prot)<-c("pathway","score","direction","fdr", "name")

sig_rna<-cbind.data.frame(rna.sub.m$pathway,-log(rna.sub.m$fdr),rna.sub.m$updown,rna.sub.m$fdr, rna.sub.m$short)
colnames(sig_rna)<-c("pathway","score","direction","fdr", "name")

sig_phos<-cbind.data.frame(phos.sub.m$pathway,-log(phos.sub.m$fdr),phos.sub.m$updown,phos.sub.m$fdr, phos.sub.m$short)
colnames(sig_phos)<-c("pathway","score","direction","fdr", "name")

sig_rna.ord<-sig_rna[order(match(sig_rna$pathway,sig_prot$pathway)),]
sig_phos.ord<-sig_phos[order(match(sig_phos$pathway,sig_prot$pathway)),]
sig_prot.ord<-sig_prot

sig_all<-rbind.data.frame(sig_prot.ord,sig_rna.ord,sig_phos.ord)
identical(as.character(as.factor(sig_prot.ord$pathway)),as.character(as.factor(sig_rna.ord$pathway)))

dtf11<-sig_all
dtf11$data<-c(rep("Global_Protein",27),rep("RNA",27),rep("Phospho",20))

dtf11$pathway<- toupper(dtf11$pathway)
dtf11$pathway<- gsub(" ", "_", dtf11$pathway)

dtf11$pathway<-gsub("KEGG_", "", dtf11$pathway)
dtf11$pathway<-gsub("DDR_", "", dtf11$pathway)
dtf11$pathway<-gsub("HALLMARK_", "", dtf11$pathway)
dtf11$pathway<-gsub("REACTOME_", "", dtf11$pathway)
dtf11$pathway<-gsub("DQ_", "", dtf11$pathway)
dtf11$pathway<-gsub("TGFb_", "", dtf11$pathway)

dtf11$pathway<-as.factor(dtf11$pathway)

dtf11$direction<-ifelse(dtf11$direction=="up",1,-1)
dtf11$sign<-dtf11$score*dtf11$direction
dtf11$direction<-as.factor(dtf11$direction)

s<-sort(dtf11$sign[dtf11$data=="Global_Protein"],ind=T,decreasing = T)
#paths<-dtf11$pathway[s$ix]
paths<- dtf11$name[s$ix]
#dtf11$pathway<-factor(dtf11$pathway,levels=paths[length(paths):1])
dtf11$name<-factor(dtf11$name,levels=paths[length(paths):1])

library(tidyverse)
library(RColorBrewer)

dtf.all<-dtf11

cols <- c(brewer.pal(n=9,name = 'Reds')[c(8,4)],brewer.pal(n=9,name = 'Greens')[c(8,4)],brewer.pal(n=9,name = 'Purples')[c(8,4)])


cols_fn <- function(data,significant)
{
  if (data == "Global_Protein" & significant == "positive") col <- cols[1]
  if (data == "Global_Protein" & significant == "negative") col <- cols[2]
  if (data == "RNA" & significant == "positive") col <- cols[3]
  if (data == "RNA" & significant == "negative") col <- cols[4]
  if (data == "Ubiquitin" & significant == "positive") col <- cols[5]
  if (data == "Ubiquitin" & significant == "negative") col <- cols[6]
  return(col)
}
dtf.all <- as_tibble(dtf.all)
dtf.all <- dtf.all %>%
  mutate(group=ifelse(fdr < 0.1, paste(data,"FDR<0.1"), paste(data,"FDR>0.1")))

p<-ggplot(dtf.all, aes(x=name, y=sign,fill=group)) +
  geom_dotplot(binaxis='y', binwidth = 1,stackratio=1.5, dotsize=5) + coord_flip() + theme_classic() + geom_hline(yintercept = 0, linetype = 2) + labs(x="Pathway",y="-log(FDR)")+
  theme(axis.text.x = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),  legend.text=element_text(colour="black",size=18), legend.title=element_text(colour="black",size=20), plot.margin=unit(c(1,-0.5,1,1), "cm")) 
dot<-p + theme(panel.background = element_rect(fill = NA),panel.grid.major.y = element_line(colour = "grey50",linetype = "dashed"),panel.grid.minor.y = element_line(colour = "grey50",linetype = "dashed")) + scale_fill_manual(values=cols) + theme(legend.position='none')


######## barplot of overlapping genes with literature ############

path<- read.csv("pathway_marcin.csv", header = T)
path1<- path[!is.na(path$short),]

path$pathway<-gsub("KEGG_", "", path$pathway)
path$pathway<-gsub("DDR_", "", path$pathway)
path$pathway<-gsub("HALLMARK_", "", path$pathway)
path$pathway<-gsub("REACTOME_", "", path$pathway)
path$pathway<-gsub("DQ_", "", path$pathway)
path$pathway<-gsub("TGF_", "", path$pathway)

path$pathway<- toupper(path$pathway)
path$pathway<- gsub(" ", "_", path$pathway)

path1$short<- paste0(path1$short," (",path1$Total.genes,")")
colnames(path1)[2]<- "name"

ord<- order(match(path1$name, levels(dtf.all$name)))

path.ord<- path1[ord,]
colnames(path.ord)[4]<- "overlap"
path.ord$name<- factor(path.ord$name, levels=(path.ord$name))

ggplot(path.ord, aes(x=name, y=prop)) +
  geom_bar(stat="identity", fill="grey") + coord_flip() + theme_classic() + labs(x="",y="Percentage")+
  theme(axis.text.x = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_blank(),axis.line.y=element_blank(), axis.ticks.y = element_blank(), 
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"), plot.margin=unit(c(1,1,1,-0.5), "cm")) + theme(legend.position = "none") 

############## pathway validation heatmap ############

path<- read.csv("pathway_valid_heat_marcin.csv")
path1<- path[,-1]
path1<- as.matrix(path1)

palette = colorRampPalette(c("grey", "red", "green","purple")) (9)
palette = colorRampPalette(c("grey", "#FF6666")) (20)
col_fun = colorRamp2(c(0,1,2,3), c("grey","red3","green4", "purple"))

Heatmap(path1, col=col_fun, rect_gp = gpar(col = "white", lwd = 2), show_column_dend = F, show_row_dend = FALSE,  cluster_rows = F, column_order = colnames(path1),row_names_side = "right",row_names_gp = gpar(fontsize = 9.5), show_row_names = F, show_column_names = F, show_heatmap_legend=F)
