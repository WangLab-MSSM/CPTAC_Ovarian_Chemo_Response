############## mappng colnames of the data to sample.anno ##############

col.name<- function(data, sample.anno){
  
  br<- which(grepl("bridge",sample.anno$Process.ID)==T)
  sample.anno<- sample.anno[-br,]

    sample.anno.id<- substr(sample.anno$Data.label..Data.set.code._.TMT.plex._.tumor.code., 11, nchar(as.character(as.factor(sample.anno$Data.label..Data.set.code._.TMT.plex._.tumor.code.)))-4)
    
  sample.anno.id<- str_replace(sample.anno.id, "-", ".")
  
  ind.br<- which(grepl("bridge", names(data))==T)
  
  data.id<- substr(names(data)[-ind.br], 11, nchar(names(data)[-ind.br])-4)
  
  if(identical(sample.anno.id, data.id) == T){
    names(data)[-ind.br]<- as.character(as.factor(sample.anno$Data.Label.2..data.set._.TMTplex._.sample.source._.patient.ID._sample.ID.))
    }else{
      print("order do not match")
  }

  return(data)
}

############## log2.ratio function ##############

log2.ratio<- function(data, ref.bridge){
  
  data[data==0]<- NA
  
  ind.ref.br<- which(grepl(ref.bridge, colnames(data))==T)
  
  data.log2=log(data,2)
  
  data.log2<- as.matrix(data.log2)
  data.log2.ratio<- as.matrix(data.log2)
  
  num.samp<- dim(data)[2]/num.tmt
  
  for(i in 1:num.tmt)
  {
    data.log2.ratio[,c((num.samp*(i-1)+1):(num.samp*i))]<- data.log2[,c((num.samp*(i-1)+1):(num.samp*i))] - data.log2[,ind.ref.br[i]]
  }
  
  data.log2.ratio.refrem<- data.log2.ratio[,-ind.ref.br]
  
  return(data.log2.ratio.refrem)
}

############## normalization function ##############
my.gnf<-function(x.m, other.bridge)
{
  if(other.bridge==T){
    ind.br<-which(grepl("bridge",colnames(x.m))==T)
    
    x.m1<- x.m[,-ind.br]
    
    m.v=apply(x.m1, 2, function(x) median(x,na.rm=T))
    
    cm=mean(m.v,na.rm=T)
    
    csd=mean(apply(abs(x.m1-matrix(m.v, nrow=nrow(x.m1), ncol=ncol(x.m1), byrow=T)), 2, function(x) median(x,na.rm=T)),na.rm=T)
  }else{
    
    m.v=apply(x.m, 2, function(x) median(x,na.rm=T))
    
    cm=mean(m.v,na.rm=T)
    
    csd=mean(apply(abs(x.m-matrix(m.v, nrow=nrow(x.m), ncol=ncol(x.m), byrow=T)), 2, function(x) median(x,na.rm=T)),na.rm=T)
  }
  
  result=apply(x.m, 2, FUN=my.gn, cm=cm, csd=csd)
}

my.gn<-function(x.v, cm, csd)
{
  temp=(x.v-median(x.v,na.rm=T))
  mad=median(abs(temp),na.rm=T)
  result=temp/mad*csd+cm
}

normalize.func<- function(data, other.bridge){
  ############ global normalization ################
  
  data.LR<- data
  
  if(other.bridge==T){
    
    data.intn=my.gnf(data.LR, other.bridge = T)  
    
  }else{
    
    data.intn=my.gnf(data.LR, other.bridge = F)  
  }
  
  return(data.intn)
}

##################### TMT outlier datapoints removal ##############
  
TMT.out.rem<- function(data, thresh = 10^(-10)){

  ind.br<- which(grepl("bridge", colnames(data))==T)
  
  data.intn.wobr<- data[,-ind.br]
  
  pval.mat<-matrix(NA,nrow(data.intn.wobr),num.tmt)
  
  for(k in 1:nrow(data.intn.wobr))
  {
    for(i in 1:num.tmt)
    {
      if(sum(is.na(data.intn.wobr[k,c((num.samples.per.tmt*(i-1)+1):(num.samples.per.tmt*i))]))==num.samples.per.tmt)
      {
        pval.mat[k,i]<-NA
      }else if(sum(is.na(data.intn.wobr[k,-c((num.samples.per.tmt*(i-1)+1):(num.samples.per.tmt*i))]))==((num.tmt-1)*(num.samples.per.tmt)))
      {
        pval.mat[k,i]<-NA
      }
      else{
        tt<-t.test(data.intn.wobr[k,c((num.samples.per.tmt*(i-1)+1):(num.samples.per.tmt*i))],data.intn.wobr[k,-c((num.samples.per.tmt*(i-1)+1):(num.samples.per.tmt*i))],var.equal=TRUE) 
        pval.mat[k,i]<-tt$p.value
      }
    }
  }
  
  m<-which(pval.mat< thresh,arr.ind = T)
  
  data.no.out<-data.intn.wobr
  for(i in 1:dim(m)[1])
  {
    j=m[i,2]
    data.no.out[m[i,1],c((num.samples.per.tmt*(j-1)+1):(num.samples.per.tmt*j))]<-NA
  }
  
  return(data.no.out)
}  

############ missing rate based filtering for imputation ################
  
missing.filt<- function(data, other.cols, sample.anno){
  
  br<- which(grepl("bridge",sample.anno$Process.ID)==T)
  sample.anno<- sample.anno[-br,]
  
  data.no<- data
  
  sen.sample<- data.no[,sample.anno$Tumor.Response=="sensitive"]
  ref.sample<- data.no[,sample.anno$Tumor.Response=="refractory"]
  
  sen.na<-apply(sen.sample,1,function(x) mean(is.na(x)))
  ref.na<-apply(ref.sample,1,function(x) mean(is.na(x)))
  
  ind.sen<-which(sen.na<0.5)
  ind.ref<-which(ref.na<0.5)
  
  ind<-union(ind.sen,ind.ref)
  
  data.filt<-data.no[ind,]
  
  d.m.des<- other.cols
  
  d.m.des.filt<-d.m.des[ind,]
  
  return(list("data"=data.filt,  "other"=d.m.des.filt))
}
  
################# check if batch correction needed
batch.check.pc<- function(data){
  
  ################# check if batch correction needed ################
  data.imp<- data
  
  data.knn<- impute.knn(as.matrix(data.imp) ,k = 10, rowmax = 0.9, colmax = 0.9, maxp = nrow(data.imp), rng.seed=362436069)
  
  data.knn1<-data.knn$data
  
  forpca<-data.knn1
  
  nm<-rep(NA,num.tmt*num.samples.per.tmt)
  for(i in 1:num.tmt)
  {
    nm[(num.samples.per.tmt*(i-1)+1):(num.samples.per.tmt*i)]<-paste0("b_",i)
  }
  
  colnames(forpca)<-make.names(nm, unique = T)
  
  colnames(forpca)<- factor(colnames(forpca), levels=colnames(forpca))
  
  fp<-t(forpca)
  cellglob.pca <- prcomp(fp)
  
  dat<-cellglob.pca$x
  dat<-data.frame(dat)
  
  dat$samples<- nm
  dat$samples<- factor(dat$samples, levels = unique(dat$samples))
  
  pc12<- ggplot(dat, aes(PC1, PC2, label = samples)) + geom_point(aes(col=samples)) + xlab("PC1") + ylab("PC2") + theme_bw() + geom_text(aes(size=3,col=samples)) +
    theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) + theme(legend.position = "none") + stat_ellipse(geom = "polygon", aes(fill=samples), alpha=0.1)
  
  pc13<- ggplot(dat, aes(PC1, PC3, label = samples)) + geom_point(aes(col=samples)) + xlab("PC1") + ylab("PC3") + theme_bw() + geom_text(aes(size=3,col=samples)) +
    theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) + theme(legend.position = "none") + stat_ellipse(geom = "polygon", aes(fill=samples), alpha=0.1)
  
  pc23<- ggplot(dat, aes(PC2, PC3, label = samples)) + geom_point(aes(col=samples)) + xlab("PC2") + ylab("PC3") + theme_bw() + geom_text(aes(size=3,col=samples)) +
    theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) + theme(legend.position = "none") + stat_ellipse(geom = "polygon", aes(fill=samples), alpha=0.1)
  
  pc14<- ggplot(dat, aes(PC1, PC4, label = samples)) + geom_point(aes(col=samples)) + xlab("PC1") + ylab("PC4") + theme_bw() + geom_text(aes(size=3,col=samples)) +
    theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) + theme(legend.position = "none") + stat_ellipse(geom = "polygon", aes(fill=samples), alpha=0.1)
  
  pc24<- ggplot(dat, aes(PC2, PC4, label = samples)) + geom_point(aes(col=samples)) + xlab("PC2") + ylab("PC4") + theme_bw() + geom_text(aes(size=3,col=samples)) +
    theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) + theme(legend.position = "none") + stat_ellipse(geom = "polygon", aes(fill=samples), alpha=0.1)
  
  pc34<- ggplot(dat, aes(PC3, PC4, label = samples)) + geom_point(aes(col=samples)) + xlab("PC3") + ylab("PC4") + theme_bw() + geom_text(aes(size=3,col=samples)) +
    theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) + theme(legend.position = "none") + stat_ellipse(geom = "polygon", aes(fill=samples), alpha=0.1)
  
 p<- ggarrange(pc12, pc13, pc14, pc23, pc24, pc34, nrow = 2, ncol = 3)
 return(p)
 
}

batch.check.icc<- function(data){

  data.imp<- data
  data.knn<- impute.knn(as.matrix(data.imp) ,k = 10, rowmax = 0.9, colmax = 0.9, maxp = nrow(data.imp), rng.seed=362436069)
  
  data.knn1<-data.knn$data
  
  icc.vec<- NULL

  for(i in 1:dim(data.knn1)[1])
  {
    data<-cbind.data.frame("x"=rep(1:num.tmt,each=num.samples.per.tmt),"y"=t(data.imp[i,]))
    colnames(data)<-c("batch","prot")
    rownames(data)<- NULL

    data$batch<- as.factor(data$batch)
    
    suppressMessages(temp<- try(ICCest(x=batch,y=prot,data=data), silent = T))

    if(class(temp)!="try-error"){
      icc.vec<-c(icc.vec,temp$ICC)
    }
    sink()
  }
  
  return(icc.vec)
}

########### Batch correction ############

get.na.ind<- function(data){
  
  data.imp<- data
  
  data.imp.vec<-c(as.matrix(data.imp))
  
  na.ind<-which(is.na(data.imp.vec))
  
  return(na.ind)
}

batch.corr<- function(data){
  
  data.imp<- data
  
  data.knn<- impute.knn(as.matrix(data.imp) ,k = 10, rowmax = 0.9, colmax = 0.9, maxp = nrow(data.imp), rng.seed=362436069)
  
  data.b<-data.knn$data
  
  na.ind<- get.na.ind(data.imp)
  
  #sink("NULL")
  
  data.nocb=ComBat(as.matrix(data.b), batch=rep(1:num.tmt,each=num.samples.per.tmt))
  
  #sink()
  
    data.nocb.vec<-c(as.matrix(data.nocb))
  data.nocb.vec[na.ind]<- NA
  
  data.nocb.out<-matrix(data.nocb.vec,nrow(data.nocb),ncol(data.nocb))
  
  colnames(data.nocb.out)<- colnames(data.imp)
  
  return(data.nocb.out)
}

append<- function(data, other.cols)
{
  data<- cbind.data.frame(other.cols, data)
  return(data)
}

########### collapse replicates ###################
# data.nocb.coll<- result$Ensemble
# 
# sample.anno<- read.csv("samples_FFPE_Aug19.csv")
# sample.anno$Tumor.Response[sample.anno$Tumor.Response=="Sensitive"]<- "sensitive"
# sample.anno$Tumor.Response[sample.anno$Tumor.Response=="Refractory"]<- "refractory"
# br<- which(grepl("bridge",sample.anno$Process.ID)==T)
# sample.anno<- sample.anno[-br,]
# 
# id<- c("D8","E5","H3","M8","R8")
# 
# rep1<-which(grepl(id[1],sample.anno$Process.ID)==T) # "p2_D8_R_Prime_OV"  "p10_D8_R_Prime_OV" "p18_D8_R_Prime_OV"
# rep2<-which(grepl(id[2],sample.anno$Process.ID)==T) # "p8_E5_R_Prime_PT"  "p16_E5_R_Prime_PT"
# rep3<-which(grepl(id[3],sample.anno$Process.ID)==T) # "p4_H3_R_Met_OM"  "p12_H3_R_Met_OM" "p19_H3_R_Met_OM"
# rep4<-which(grepl(id[4],sample.anno$Process.ID)==T) # "p6_M8_S_Prime_other"  "p14_M8_S_Prime_other"
# rep5<-which(grepl(id[5],sample.anno$Process.ID)==T) # "p1_R8_R_Prime_PT"  "p9_R8_R_Prime_PT"  "p17_R8_R_Prime_PT"
# 
# avg.rep1<-apply(data.nocb.coll[,rep1],1,function(x) mean(x,na.rm=T))
# avg.rep2<-apply(data.nocb.coll[,rep2],1,function(x) mean(x,na.rm=T))
# avg.rep3<-apply(data.nocb.coll[,rep3],1,function(x) mean(x,na.rm=T))
# avg.rep4<-apply(data.nocb.coll[,rep4],1,function(x) mean(x,na.rm=T))
# avg.rep5<-apply(data.nocb.coll[,rep5],1,function(x) mean(x,na.rm=T))
# 
# data.nocb.coll[,rep1[1]]<-avg.rep1
# data.nocb.coll[,rep2[1]]<-avg.rep2
# data.nocb.coll[,rep3[1]]<-avg.rep3
# data.nocb.coll[,rep4[1]]<-avg.rep4
# data.nocb.coll[,rep5[1]]<-avg.rep5
# 
# data.nocb1.coll<-data.nocb.coll[,-c(rep1[-1],rep2[-1],rep3[-1],rep4[-1],rep5[-1])]
# 
# rownames(data.nocb1.coll)<-batch.corrected.missing.filt.share$Index
# 
# save(data.nocb1.coll, file="FFPE_imp_coll.RData")
# 
# 
# 
