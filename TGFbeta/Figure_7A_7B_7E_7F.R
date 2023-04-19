t<- read.csv("FFPE_protein_TGFb_ALTEJ_scores.tsv",sep="\t")
t<- read.csv("FFPE_RNA_TGFb_ALTEJ_scores.tsv",sep="\t")
t<- read.csv("Frozen_protein_TGFb_ALTEJ_scores.tsv",sep="\t")
t<- read.csv("FFPE_validation_TGFb_ALTEJ_scores.tsv",sep="\t")
t<- read.csv("CPTAC2_protein_TGFb_ALTEJ_scores.tsv",sep="\t")

ggplot(t, aes(TGFb,Alt_EJ)) + geom_point() + xlab("TGFb score") + ylab("Alt-EJ score") + geom_smooth(method='lm', formula= Alt_EJ~TGFb) + theme_bw() + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 

cor.test(t$TGFb, t$Alt_EJ)

############## 7B ##############
sample.anno<- read.csv("samples_FFPE_Aug19.csv")
br<- which(grepl("bridge",sample.anno$Process.ID)==T)
sample.anno<- sample.anno[-br,]

sample.anno$Tumor.Response[sample.anno$Tumor.Response=="Sensitive"]<- "sensitive"
sample.anno$Tumor.Response[sample.anno$Tumor.Response=="Refractory"]<- "refractory"


id<-c("D8","E5","H3","M8","R8")

rep1<-which(sample.anno$Process.ID=="D8") 
rep2<-which(sample.anno$Process.ID=="E5") 
rep3<-which(sample.anno$Process.ID=="H3") 
rep4<-which(sample.anno$Process.ID=="M8") 
rep5<-which(sample.anno$Process.ID=="R8") 

sample.anno.worep<-sample.anno[-c(rep1[-1],rep2[-1],rep3[-1],rep4[-1],rep5[-1]),]

sen<-which(sample.anno.worep$Tumor.Response=="sensitive")

ref<-which(sample.anno.worep$Tumor.Response=="refractory")

senref<-rep(NA,158)

senref[sen]<-rep(1,length(sen)) ## 91
senref[ref]<-rep(0,length(ref)) ## 67
table(senref)

cov<- read.csv("covariates.csv", header = T)

cov.ord<- cov[order(match(cov$Process.ID,sample.anno.worep$Process.ID)),]

ov<-which(cov.ord$Tumor.Location.Group=="OV")
om<-which(cov.ord$Tumor.Location.Group=="OM")

location<-rep("A",158)
location[ov]<-rep("ov",length(ov)) 
location[om]<-rep("om",length(om)) 

fh<- which(sample.anno.worep$Sample.Source=="FH")
mc<- which(sample.anno.worep$Sample.Source=="MC")

site<- rep("A",158)
site[fh]<- rep("FH",length(fh))
site[mc]<- rep("MC",length(mc))

age<- cov.ord$Patient.Age
neo<- cov.ord$Neo.adjuvant

load("tumor_purity.RData")
tp<- as.numeric(as.matrix(tp1[,2]))

t<- read.csv("prot_assoc/FFPE_protein_TGFb_ALTEJ_scores.tsv",sep="\t")

df<- cbind.data.frame(samples = sample.anno.worep$Data.Label.2..data.set._.TMTplex._.sample.source._.patient.ID._sample.ID., t,"response"=senref, "loc"=location,"site"=site,"age"=age,"neo"=neo, "tp"=tp)

class<- read.csv("km7prot/km7prot.k=5.consensusClass.csv", header = F)
colnames(class)[1]<- "samples"
colnames(class)[2]<- "clusters"

df<- df[,-1]
m<- merge(df,class,by="samples")
m$clusters<- as.factor(m$clusters)
m$response<- as.factor(m$response)          

mod1<- lm(beta_alt ~ clusters + loc + site + age + neo + tp, data=m) ## change beta_alt to TGFb and Alt_EJ
mod2<- lm(beta_alt ~ loc + site + age + neo + tp, data=m) ## change beta_alt to TGFb and Alt_EJ

anova(mod1,mod2, test = "Chisq") 

ggplot(m, aes(x=clusters, y=beta_alt, fill=clusters)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = clusters), size = 2, shape = 21, position = position_jitterdodge()) + theme_bw() + ylab("Beta_alt score") + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") ### change beta_alt to TGFb and Alt_EJ

############## 7F #############
m.sub<- m[m$clusters==4,]

ggplot(m.sub, aes(x=clusters, y=TGFb, fill=response)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = response), size = 2, shape = 21, position = position_jitterdodge()) + theme_bw() + ylab("TGFb score") + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 

t.test(m.sub$TGFb[m$response==1], m.sub$TGFb[m$response==0])

path<- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
t<- hall[which(hall$pathway %in% path),]
t<- t[,-1]
pn<- t$pathway
t<- as.matrix(t)
rownames(t)<- path

df<- cbind.data.frame(samples = sample.anno.worep$Data.Label.2..data.set._.TMTplex._.sample.source._.patient.ID._sample.ID., t,"response"=senref, "loc"=location,"site"=site,"age"=age,"neo"=neo, "tp"=tp)

df<- df[,-1]
m<- merge(df,class,by="samples")
m$clusters<- as.factor(m$clusters)
m$response<- as.factor(m$response)          
m.sub<- m[m$clusters==4,]

ggplot(m.sub, aes(x=clusters, y=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, fill=response)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = response), size = 2, shape = 21, position = position_jitterdodge()) + theme_bw() + ylab("EMT score") + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 

############## 7E #############
t<- read.csv("FFPE_RNA_TGFb_ALTEJ_scores.tsv",sep="\t")

sample.anno<-read.csv("rna_assoc/samples_RNA.csv", header = T)
sample.anno$Tumor.response[sample.anno$Tumor.response=="Sensitive"]<- "sensitive"
sample.anno$Tumor.response[sample.anno$Tumor.response=="Refractory"]<- "refractory"
sample.anno$RNA.ID<- as.character(sample.anno$RNA.ID)

sample.anno$RNA.ID[c(7,18,19,21,27)]<- c("R69", "R80", "R82", "R85", "R93")

sample.anno.sub<- sample.anno[sample.anno$RNA.ID %in% colnames(rna),]

#sample.anno.sub<- sample.anno[sample.anno$RNA.ID %in% ba$samples,]

sample.anno.sub.ord<- sample.anno.sub[order(match(sample.anno.sub$RNA.ID, colnames(rna))),]

#sample.anno.sub.ord<- sample.anno.sub[order(match(sample.anno.sub$RNA.ID, ba$samples)),]

identical(sample.anno.sub.ord$RNA.ID, colnames(rna)) # TRUE

sen<-which(sample.anno.sub.ord$Tumor.response=="sensitive")

ref<-which(sample.anno.sub.ord$Tumor.response=="refractory")

senref<-rep(NA,106)

senref[sen]<-rep(1,length(sen)) ## 58
senref[ref]<-rep(0,length(ref)) ## 48
table(senref)

neo<- sample.anno.sub.ord$Neo.adjuvant
neo1<- ifelse(neo=="yes",1,0)

age<- sample.anno.sub.ord$Patient.Age
om<-which(sample.anno.sub.ord$Tumor.Location.Group=="OM")

ov<-which(sample.anno.sub.ord$Tumor.Location.Group=="OV")

location1<-rep("A",106)
location1[om]<-rep("OM",length(om)) 
location1[ov]<-rep("OV",length(ov)) 

fh<- which(sample.anno.sub.ord$Sample.Source=="FHCRC")
mc<- which(sample.anno.sub.ord$Sample.Source=="Mayo")

site<- rep("A",106)
site[fh]<- rep("FH",length(fh))
site[mc]<- rep("MC",length(mc))

tp1<- read.table("rna_assoc/rna_tp.txt", header=T)
tp<- as.numeric(tp1$ESTIMATE)

df<- cbind.data.frame("prot_id" = sample.anno.worep$Data.Label.2..data.set._.TMTplex._.sample.source._.patient.ID._sample.ID., t,"response"=senref, "loc"=location,"site"=site,"age"=age,"neo"=neo1, "tp"=tp)

load("prot_rna_map_new.RData")
m.ord$samples<- as.character(as.factor(m.ord$samples))
m.ord$samples[62]<- "R22"

df1<- merge(m.ord, df, by="prot_id")

class<- read.csv("km7prot/km7prot.k=5.consensusClass.csv", header = F)
colnames(class)[1]<- "samples"
colnames(class)[2]<- "clusters"

m<- merge(df1,class,by="samples")
m$clusters<- as.factor(m$clusters)
m$response<- as.factor(m$response)          

m.sub<- m[m$clusters==4,]

ggplot(m.sub, aes(x=clusters, y=TGFb, fill=response)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = response), size = 2, shape = 21, position = position_jitterdodge()) + theme_bw() + ylab("TGFb score") + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 
