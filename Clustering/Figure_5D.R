anno<- read.table("mrm_sample_tmt.txt", header = T, sep="\t")

mrm<- read.table("mrm_data.txt", header=T, sep="\t")
mrm1<- mrm[,-c(1:5)]
other<- mrm[,c(1:5)]

row.na<- apply(mrm1, 1, function(x) mean(is.na(x)))

mrm3<- mrm1[which(row.na < 0.5),]
other3<- other[which(row.na < 0.5),]

row.na.filt<- apply(mrm3, 1, function(x) mean(is.na(x)))

hist(row.na, xlab="missing rate per peptide")

anno$MRM_ID<- str_replace(anno$MRM_ID, " ", ".")

identical(anno$MRM_ID, colnames(mrm3)) # TRUE

mrm4<- t(apply(mrm3, 1, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T)))

ung<- unique(other3$Gene)

pval<- NULL; fc<- NULL 
prot<- NULL
for(i in 1:length(unique(other3$Gene))){
  foo<- mrm4[which(other3$Gene==ung[i]),]
  
  if(is.null(dim(foo))){
    foo.avg<- as.numeric(foo)
  }else{
    foo.avg<- as.numeric(apply(foo, 2, function(x) mean(x, na.rm=T)))
  }
  prot<- rbind(prot, foo.avg)
  cl3<- foo.avg[which(anno$Cluster==3)]
  rest<- foo.avg[which(anno$Cluster!=3)]
  
  temp<- try(t.test(cl3, rest), silent = T)
  
  if(class(temp)!="try-error"){
    pval<- c(pval, temp$p.value)
    fc<- c(fc, 2^(temp$estimate[1] - temp$estimate[2]))
  }else{
    pval<- c(pval, NA)
    fc<- c(fc, NA)
  }
}

rownames(prot)<- ung
tab<- cbind.data.frame("Gene"=ung, "pval"=pval, "fc"=fc)

########### boxplot of proteins based on FFPE TMT data ##########
load("FFPE_imputed_coll.RData")
data<-data.nocb1.coll
other<- data[,c(1:9)]
data1<- data[,-c(1:9)]
data1<- as.matrix(data1)
rownames(data1)<- data$symbol

data2<- data1[,which(colnames(data1) %in% anno$TMT_ID)]
data3<- data2[,order(match(colnames(data2), anno$TMT_ID))]

data4<- data3[which(rownames(data3) %in% unique(other$Gene)),]

identical(colnames(data3),anno$TMT_ID)#TRUE 

pval<- NULL; fc<- NULL 
for(i in 1:nrow(data4)){
  cl3<- as.numeric(data4[i,which(anno$Cluster==3)])
  rest<- as.numeric(data4[i,which(anno$Cluster!=3)])
  
  temp<- try(t.test(cl3, rest), silent = T)
  
  if(class(temp)!="try-error"){
    pval<- c(pval, temp$p.value)
    fc<- c(fc, 2^(temp$estimate[1] - temp$estimate[2]))
  }else{
    pval<- c(pval, NA)
    fc<- c(fc, NA)
  }
}

tab<- cbind.data.frame("Gene"=rownames(data4), "pval"=pval, "fc"=fc)

anno$group<- ifelse(anno$Cluster==3, 1, 0)

data.gene<- data3[which(rownames(data3) %in% c("AIFM1","COX6C", "COX7A2", "GLUD1", "ATP5F1B", "ATP5F1C", "HADH", "CPOX")),]
data.gene.sd<- t(apply(data.gene,1,function(x) (x-mean(x))/sd(x)))

df<- cbind.data.frame("samples"=colnames(data.gene.sd),t(data.gene.sd))
rownames(df)<- NULL

colnames(anno)[2]<- "samples"

class<- read.csv("km7prot/km7prot.k=5.consensusClass.csv", header = F)
colnames(class)[1]<- "samples"
m<- merge(df,class,by="samples")
m<- merge(df,anno,by="samples")

colnames(m)[14]<- "clusters"
m$clusters<- as.factor(m$clusters)

dat<- cbind.data.frame("prot"=c(m[,2], m[,3], m[,4], m[,5], m[,6], m[,7], m[,8], m[,9]), "genes"=rep(colnames(m)[2:9], each=102), "clusters"=rep(m$clusters,8))

ggplot(dat, aes(x=genes, y=prot, fill=clusters)) + geom_boxplot() + theme_bw() + ylab("Protein abundance") + theme(axis.text.x = element_text(colour="black",size=20,angle=45,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"))


########### boxplot of proteins based on MRM data ############
df.prot<- prot[which(rownames(prot) %in% c("AIFM1","COX6C", "COX7A2", "GLUD1", "ATP5F1B", "ATP5F1C", "HADH", "CPOX")),]
df.prot.sd<- t(apply(df.prot, 1, function(x) (x-mean(x))/sd(x)))

dat1<- c(df.prot.sd[1,], df.prot.sd[2,], df.prot.sd[3,], df.prot.sd[4,], df.prot.sd[5,],df.prot.sd[6,],df.prot.sd[7,],df.prot.sd[8,])
dat<- cbind.data.frame("mrm"=dat1,"genes"=rep(rownames(df.prot), each=102), "clusters" = rep(anno$Cluster, 8))
dat$clusters<- as.factor(dat$clusters)

ggplot(dat, aes(x=genes, y=mrm, fill=clusters)) + geom_boxplot() + theme_bw() + ylab("MRM data") + ylim(c(-2,3)) + theme(axis.text.x = element_text(colour="black",size=20,angle=45,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"))

########### boxplot of proteins based on CPTAC2 data ##########

data<- CPTAC.Test.Final

data.gene<- data1[which(rownames(data1) %in% c("AIFM1","COX6C", "COX7A2", "GLUD1", "ATP5F1B", "ATP5F1C", "HADH", "CPOX")),]

data.gene.sd<- t(apply(data.gene,1,function(x) (x-mean(x))/sd(x)))

df<- cbind.data.frame("samples"=colnames(data.gene.sd),t(data.gene.sd))
rownames(df)<- NULL

class<- read.csv("CPTAC2/CPTAC_PAM.txt", sep="\t",header = T)

class<- class[1:174,]
colnames(class)[1]<- "samples"
m<- merge(df,class,by="samples")
colnames(m)[10]<- "clusters"
m$clusters<- as.factor(m$clusters)

dat<- cbind.data.frame("prot"=c(m[,2], m[,3], m[,4], m[,5], m[,6], m[,7], m[,8], m[,9]), "genes"=rep(colnames(m)[2:9], each=174), "clusters"=rep(m$clusters,8))

ggplot(dat, aes(x=genes, y=prot, fill=clusters)) + geom_boxplot() + theme_bw() + ylab("Protein abundance") + theme(axis.text.x = element_text(colour="black",size=20,angle=45,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"))

pval<- NULL; fc<- NULL 
for(i in 2:9){
  cl3<- as.numeric(m[which(m$clusters==3),i])
  rest<- as.numeric(m[which(m$clusters!=3),i])
  
  temp<- t.test(cl3, rest)
  
  pval<- c(pval, temp$p.value)
  fc<- c(fc, 2^(temp$estimate[1] - temp$estimate[2]))
}

tab<- cbind.data.frame("Gene"=colnames(m)[2:9], "pval"=pval, "fc"=fc)

write.csv(tab, file = "TMT158_cl3_vs_rest.csv")

########### roc cl3 vs rest ###########
df.prot<- prot[which(rownames(prot) %in% c("AIFM1","COX6C", "COX7A2", "GLUD1", "ATP5F1B", "ATP5F1C", "HADH", "CPOX")),]
df.prot.sd<- t(apply(df.prot, 1, function(x) (x-mean(x))/sd(x)))

avg.score<- apply(df.prot.sd, 2, function(x) mean(x, na.rm=T))
anno$group<- ifelse(anno$Cluster==3, 1, 0)

tab<- cbind.data.frame("samples"=colnames(mrm3), "pred.scores"=avg.score)
write.csv(tab, file = "pred_scores_mrm_cluster3.csv")

roc<- roc(anno$group, avg.score, ci = TRUE, direction="<")

ciobj1<-ci.se(roc, specificities = seq(0,1,.01))
dat.ci1 <- data.frame(x = as.numeric(rownames(ciobj1)),
                      lower = ciobj1[, 1],
                      upper = ciobj1[, 3])

gtrain<-roc
ggroc(gtrain,alpha = 0.5, colour = "green4", linetype = 1, size = 2) +  theme_bw() + annotate(geom="text", x=0.2, y=0.4, label="AUC = 0.84", color="red",size=7) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))  + geom_ribbon(data = dat.ci1, aes(x = x, ymin = lower, ymax = upper), fill = "green", alpha= 0.2)  
