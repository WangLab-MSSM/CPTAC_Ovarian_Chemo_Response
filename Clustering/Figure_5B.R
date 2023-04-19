mat<- read.table("data_alignment/CPTAC2_FFPE_align/pathway_scores_retro/CPTAC2_retro_pathway.txt", header = T)

mat.ord2.sd<- apply(data.prot,1,function(x) (x-mean(x))/sd(x))

library(ConsensusClusterPlus)

consensus = ConsensusClusterPlus(t(mat.ord2.sd),maxK=7,reps=100,pItem=0.8,pFeature=1,title="km7retro",distance="euclidean",clusterAlg="km",plot = "pdf", writeTable = T)

class<- read.csv("km7retro/km7retro.k=5.consensusClass.csv", header = F)
colnames(class)[1]<- "samples"
colnames(class)[2]<- "clusters"

c1<- class$V1[class$V2==1]
c2<- class$V1[class$V2==2]
c3<- class$V1[class$V2==3]
c4<- class$V1[class$V2==4]
c5<- class$V1[class$V2==5]

class.ord<- c(which(colnames(t(mat.ord2.sd)) %in% c1),which(colnames(t(mat.ord2.sd)) %in% c2),which(colnames(t(mat.ord2.sd)) %in% c3),which(colnames(t(mat.ord2.sd)) %in% c4), which(colnames(t(mat.ord2.sd)) %in% c5))

mat3<- t(mat.ord2.sd)[,class.ord]

mat5<- cbind.data.frame("samples"=colnames(mat3),t(mat3))
rownames(mat5)<- NULL

m.retro<- merge(mat5,class, by="samples")
m.retro.sub<- m.retro

m.retro1<- m.retro.sub[which(m.retro.sub$clusters==1),]
m.retro2<- m.retro.sub[which(m.retro.sub$clusters==2),]
m.retro3<- m.retro.sub[which(m.retro.sub$clusters==3),]
m.retro4<- m.retro.sub[which(m.retro.sub$clusters==4),]
m.retro5<- m.retro.sub[which(m.retro.sub$clusters==5),]

retro.c1.mean<- apply(m.retro1[,2:151],2,mean)
retro.c2.mean<- apply(m.retro2[,2:151],2,mean)
retro.c3.mean<- apply(m.retro3[,2:151],2,mean)
retro.c4.mean<- apply(m.retro4[,2:151],2,mean)
retro.c5.mean<- apply(m.retro5[,2:151],2,mean)

retro.cluster.mean<- cbind.data.frame("pathways"=rownames(mat3),"c1.mean"=retro.c1.mean,"c2.mean"=retro.c2.mean,"c3.mean"=retro.c3.mean,"c4.mean"=retro.c4.mean,"c5.mean"=retro.c5.mean)
rownames(retro.cluster.mean)<- NULL

######## FFPE Discovery cluster mean ###########
kegg<-read_tsv("prot_assoc/KEGG_ssgsea.tsv",col_names = TRUE,quote = "\"")
hall<-read_tsv("prot_assoc/Hallmark_ssgsea.tsv",col_names = TRUE,quote = "\"")
reac<-read_tsv("prot_assoc/Reactome_ssgsea.tsv",col_names = TRUE,quote = "\"")
ddr<-read_tsv("prot_assoc/DDR_ssgsea.tsv",col_names = TRUE,quote = "\"")
dq<-read_tsv("prot_assoc/DQ_ssgsea.tsv",col_names = TRUE,quote = "\"")
tgfb<-read_tsv("prot_assoc/TGFb_ssgsea.tsv",col_names = TRUE,quote = "\"")

for(i in 1:dim(kegg)[1]){
  kegg[i,1]<- as.character(kegg[i,1])
  kegg[i,1]<- paste0("KEGG_",kegg[i,1])
}
for(i in 1:dim(reac)[1]){
  reac[i,1]<- as.character(reac[i,1])
  reac[i,1]<- paste0("Reactome_",reac[i,1])
}
for(i in 1:dim(ddr)[1]){
  ddr[i,1]<- as.character(ddr[i,1])
  ddr[i,1]<- paste0("DDR_",ddr[i,1])
}
for(i in 1:dim(dq)[1]){
  dq[i,1]<- as.character(dq[i,1])
  dq[i,1]<- paste0("DQ_",dq[i,1])
}
for(i in 1:dim(tgfb)[1]){
  tgfb[i,1]<- as.character(tgfb[i,1])
  tgfb[i,1]<- paste0("TGFb_",tgfb[i,1])
}

#path<-read.csv("pathway_ssgsea.csv")

pathway<- read.csv("prot_assoc/tumor_purity/FFPE_global_pathway.tsv", sep="\t", header=T)
path<- pathway[pathway$fdr < 0.01,]

#################################  
kegg.sub<- kegg[which(kegg$pathway %in% path$pathways),]
hall.sub<- hall[which(hall$pathway %in% path$pathways),]
reac.sub<- reac[which(reac$pathway %in% path$pathways),]
ddr.sub<- ddr[which(ddr$pathway %in% path$pathways),]
dq.sub<- dq[which(dq$pathway %in% path$pathways),]
tgfb.sub<- tgfb[which(tgfb$pathway %in% path$pathways),]

mat<- rbind.data.frame(kegg.sub,hall.sub,reac.sub,ddr.sub,dq.sub,tgfb.sub)
mat1<- mat[,-1]
mat1<- as.matrix(mat1)

rownames(mat1)<- mat$pathway

mat.ord2.sd<- apply(mat1,1,function(x) (x-mean(x))/sd(x))

class<- read.csv("km7prot/km7prot.k=5.consensusClass.csv", header = F)
colnames(class)[1]<- "samples"
colnames(class)[2]<- "clusters"

c1<- class$V1[class$V2==1]
c2<- class$V1[class$V2==2]
c3<- class$V1[class$V2==3]
c4<- class$V1[class$V2==4]
c5<- class$V1[class$V2==5]

class.ord<- c(which(colnames(t(mat.ord2.sd)) %in% c1),which(colnames(t(mat.ord2.sd)) %in% c2),which(colnames(t(mat.ord2.sd)) %in% c3),which(colnames(t(mat.ord2.sd)) %in% c4), which(colnames(t(mat.ord2.sd)) %in% c5))

mat3<- t(mat.ord2.sd)[,class.ord]

mat5<- cbind.data.frame("samples"=colnames(mat3),t(mat3))
rownames(mat5)<- NULL

m.ffpe<- merge(mat5,class, by="samples")
m.ffpe.sub<- m.ffpe

m.ffpe1<- m.ffpe.sub[which(m.ffpe.sub$clusters==1),]
m.ffpe2<- m.ffpe.sub[which(m.ffpe.sub$clusters==2),]
m.ffpe3<- m.ffpe.sub[which(m.ffpe.sub$clusters==3),]
m.ffpe4<- m.ffpe.sub[which(m.ffpe.sub$clusters==4),]
m.ffpe5<- m.ffpe.sub[which(m.ffpe.sub$clusters==5),]

ffpe.c1.mean<- apply(m.ffpe1[,2:151],2,mean)
ffpe.c2.mean<- apply(m.ffpe2[,2:151],2,mean)
ffpe.c3.mean<- apply(m.ffpe3[,2:151],2,mean)
ffpe.c4.mean<- apply(m.ffpe4[,2:151],2,mean)
ffpe.c5.mean<- apply(m.ffpe5[,2:151],2,mean)

FFPE.cluster.mean<- cbind.data.frame("pathways"=rownames(mat3),"c1.mean"=ffpe.c1.mean,"c2.mean"=ffpe.c2.mean,"c3.mean"=ffpe.c3.mean,"c4.mean"=ffpe.c4.mean,"c5.mean"=ffpe.c5.mean)
rownames(ffpe.cluster.mean)<- NULL

ord<- order(match(retro.cluster.mean$pathways,FFPE.cluster.mean$pathways))

retro.cluster.mean.ord<- retro.cluster.mean[ord,]

identical(as.character(as.factor(FFPE.cluster.mean$pathways)), as.character(as.factor(retro.cluster.mean$pathways))) #TRUE

cor.test(FFPE.cluster.mean$c1.mean, retro.cluster.mean.ord$c1.mean)
cor.test(FFPE.cluster.mean$c1.mean, retro.cluster.mean.ord$c2.mean)
cor.test(FFPE.cluster.mean$c1.mean, retro.cluster.mean.ord$c3.mean)
cor.test(FFPE.cluster.mean$c1.mean, retro.cluster.mean.ord$c4.mean)
cor.test(FFPE.cluster.mean$c1.mean, retro.cluster.mean.ord$c5.mean)

cor.test(FFPE.cluster.mean$c2.mean, retro.cluster.mean.ord$c1.mean)
cor.test(FFPE.cluster.mean$c2.mean, retro.cluster.mean.ord$c2.mean)
cor.test(FFPE.cluster.mean$c2.mean, retro.cluster.mean.ord$c3.mean)
cor.test(FFPE.cluster.mean$c2.mean, retro.cluster.mean.ord$c4.mean)
cor.test(FFPE.cluster.mean$c2.mean, retro.cluster.mean.ord$c5.mean)

cor.test(FFPE.cluster.mean$c3.mean, retro.cluster.mean.ord$c1.mean)
cor.test(FFPE.cluster.mean$c3.mean, retro.cluster.mean.ord$c2.mean)
cor.test(FFPE.cluster.mean$c3.mean, retro.cluster.mean.ord$c3.mean)
cor.test(FFPE.cluster.mean$c3.mean, retro.cluster.mean.ord$c4.mean)
cor.test(FFPE.cluster.mean$c3.mean, retro.cluster.mean.ord$c5.mean)

cor.test(FFPE.cluster.mean$c4.mean, retro.cluster.mean.ord$c1.mean)
cor.test(FFPE.cluster.mean$c4.mean, retro.cluster.mean.ord$c2.mean)
cor.test(FFPE.cluster.mean$c4.mean, retro.cluster.mean.ord$c3.mean)
cor.test(FFPE.cluster.mean$c4.mean, retro.cluster.mean.ord$c4.mean)
cor.test(FFPE.cluster.mean$c4.mean, retro.cluster.mean.ord$c5.mean)

cor.test(FFPE.cluster.mean$c5.mean, retro.cluster.mean.ord$c1.mean)
cor.test(FFPE.cluster.mean$c5.mean, retro.cluster.mean.ord$c2.mean)
cor.test(FFPE.cluster.mean$c5.mean, retro.cluster.mean.ord$c3.mean)
cor.test(FFPE.cluster.mean$c5.mean, retro.cluster.mean.ord$c4.mean)
cor.test(FFPE.cluster.mean$c5.mean, retro.cluster.mean.ord$c5.mean)
