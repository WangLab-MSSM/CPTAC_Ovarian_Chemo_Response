library(ggplot2)
library(ggpubr)
prot<- read.csv("FFPE_Protein_Assoc_May18.tsv", sep="\t", header=T)
rna<- read.csv("out_assoc_RNA_site_tp_new.tsv", sep="\t", header=T)
phos<- read.csv("FFPE_Phospho_Assoc_May18.tsv", sep="\t", header=T)
colnames(rna)[1]<- "symbol"

m<- merge(prot[,c(2,10:12)],rna[,c(1,10:12)], by="symbol")
m1<- merge(m, phos[,c(2,10:12)], by="symbol")
m1<- m[complete.cases(m1),]

colnames(m1)[-1]<- c("FC_prot","pval_prot", "FDR_prot","FC_RNA","pval_RNA", "FDR_RNA","FC_phos", "pval_phos", "FDR_phos")

protp<- ifelse(m1$FC_prot>1,-log(m1$pval_prot)/log(10),log(m1$pval_prot)/log(10))
rnap<- ifelse(m1$FC_RNA>1,-log(m1$pval_RNA)/log(10),log(m1$pval_RNA)/log(10))
phosp<- ifelse(m1$FC_phos>1,-log(m1$pval_phos)/log(10),log(m1$pval_phos)/log(10))

dat<- cbind.data.frame("protp"=protp,"rnap"=rnap,"FDR_prot"=m$FDR_prot,"FDR_RNA"=m$FDR_RNA,"phosp"=phosp,"FDR_phos"=m1$FDR_phos,"genes"=m$symbol)

dat$genes.text<- rep("",dim(dat)[1])
ind.prot<- which(dat$FDR_prot < 0.1)
dat$genes.text[ind.prot]<- as.character(dat$genes[ind.prot])
ind.rna<- which(dat$FDR_RNA < 0.1)
dat$genes.text[ind.rna]<- as.character(dat$genes[ind.rna])
ind.phos<- which(dat$FDR_phos < 0.1)
dat$genes.text[ind.phos]<- as.character(dat$genes[ind.phos])

dat$col<- rep("grey", dim(dat)[1])
dat$col[c(ind.prot,ind.rna)]<- "black"
  
umut<- c("PPIL1", "MCM3", "DEK", "E2F3", "NUP153", "MCUR1", "ZSCAN16", "SRPK1", "ABCF1", "DTNBP1", "RNF144B", "KDM1B", "TPMT", "COA4", "UHRF1BP1", "TXNDC5", "TGM2")

xiaoyu<- c("LSM2", "AMY1C", "CARMIL1")

dat$symb<- rep("white",dim(dat)[1])
dat$symb[which(dat$genes %in% umut)]<- "red"
  dat$symb[which(dat$genes %in% xiaoyu)]<- "blue"
    
  dat$col[dat$genes.text!=""]<- dat$symb[dat$genes.text!=""]
  dat$col[which(dat$col=="white")]<- "black"
    
  ccdc167<- cbind.data.frame("protp"=4.341284,"rnap"=2.324206,"FDR_prot"=NA,"FDR_RNA"=NA,"genes"="CCDC167","genes.text"="CCDC167","col"="red","symb"="red")
  
  eps8l2<- cbind.data.frame("protp"=-4.494635,"rnap"=-1.569355,"FDR_prot"=NA,"FDR_RNA"=NA,"genes"="EPS8L2_s570","genes.text"="EPS8L2_s570","col"="black","symb"="white")
  
  dat1<- rbind.data.frame(dat,ccdc167,eps8l2)

  ggplot(dat1, aes(rnap, protp,col=col, label=genes.text))+ 
    geom_point(aes(col=col)) + 
    scale_color_manual(values=c("black","grey","red")) +
    geom_vline(xintercept = 0, linetype = 2)+ 
    geom_hline(yintercept = 0, linetype = 2)+
    geom_label_repel(aes(fill=symb), colour="black", segment.colour="black")+ scale_fill_manual(values=c("#FF9999","white")) +
    xlab(" signed -log_10(pvalue)(RNA)") + ylab("signed -log_10(pvalue)(Protein or Phospho)") + theme_bw() + 
    theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
          axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) + theme(legend.position = "none")
  