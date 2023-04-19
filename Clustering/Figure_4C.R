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

for(k in 1:100){
  print(k)
  
  pathway<- read.csv("FFPE_global_pathway.tsv", sep="\t", header=T)
  s<- sample(c(setdiff(1:2012,2001)), 150, replace = FALSE)
  path<- pathway[s,]
  
  kegg.sub<- kegg[which(kegg$pathway %in% path$pathway),]
  hall.sub<- hall[which(hall$pathway %in% path$pathway),]
  reac.sub<- reac[which(reac$pathway %in% path$pathway),]
  ddr.sub<- ddr[which(ddr$pathway %in% path$pathway),]
  dq.sub<- dq[which(dq$pathway %in% path$pathway),]
  tgfb.sub<- tgfb[which(tgfb$pathway %in% path$pathway),]
  
  mat<- rbind.data.frame(kegg.sub,hall.sub,reac.sub,ddr.sub,dq.sub,tgfb.sub)
  mat1<- mat[,-1]
  mat1<- as.matrix(mat1)
  
  rownames(mat1)<- mat$pathway
  
  mat2<- t(apply(mat1,1,function(x) (x-mean(x))/sd(x)))
  
  consensus = ConsensusClusterPlus(mat2,maxK=7,reps=100,pItem=0.8,pFeature=1,title=paste0("km7random",k),distance="euclidean",clusterAlg="km",plot = "pdf", writeTable = T)
  
  class<- read.csv(paste0("km7random",k,"/km7random",k,".k=5.consensusClass.csv"), header = F)
  colnames(class)[1]<- "samples"
  colnames(class)[2]<- "clusters"
  
  data<- cbind.data.frame("samples"=colnames(mat2), t(mat2))
  rownames(data)<- NULL
  
  m<- merge(data,class, by="samples")
  m$clusters<- as.factor(m$clusters)
  
  rat.vec<- NULL
  for(i in 2:151){
    fit = lm(m[,i] ~ m$clusters)
    w.var<- anova(fit)["Residuals", "Mean Sq"]
    b.var<- anova(fit)["m$clusters", "Mean Sq"]
    
    rat<- w.var/b.var
    rat.vec<- c(rat.vec,rat)
  }
  
  save(rat.vec, file = paste0("random150_invest/wbvarratio_random150/wbvarratio_random150_",k,".RData"))
}

##############################################
rat.vec.vec<- NULL
rat.mat<- NULL
for(k in 1:100){
  load(paste0("wbvarratio_random150_",k,".RData"))
  
  rat.vec.vec<- c(rat.vec.vec,rat.vec)
  rat.mat<- cbind(rat.mat,rat.vec)
}
rat.mat.vec<- apply(rat.mat,2,mean)

ggplot(tab,aes(x=rat.mat.vec, fill=group)) + geom_histogram(aes(fill=group)) + xlab("Within/Between cluster variances") + ylab("Density")  + xlim(c(0.05,.27)) + theme_bw() + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"), axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"), legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "none")