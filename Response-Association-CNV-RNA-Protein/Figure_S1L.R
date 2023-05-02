###### FFPE discovery ########
gene<- c("BCL2L1")

load("FFPE_imputed_coll.RData")

data<-data.nocb1.coll
other<- data[,c(1:9)]
data1<- data[,-c(1:9)]
data1<- as.matrix(data1)
rownames(data1)<- data$symbol

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

identical(as.character(as.factor(sample.anno.worep$Data.Label.2..data.set._.TMTplex._.sample.source._.patient.ID._sample.ID.)), colnames(data.nocb.coll.nopool)) # TRUE

sen<-which(sample.anno.worep$Tumor.Response=="sensitive")

ref<-which(sample.anno.worep$Tumor.Response=="refractory")

senref<-rep(NA,158)

senref[sen]<-rep(1,length(sen)) ## 91
senref[ref]<-rep(0,length(ref)) ## 67
table(senref)

data.prot<- data1[which(rownames(data1) %in% gene),]  
df.prot<- cbind.data.frame(t(data.prot),"Response"=senref)
df.prot$Response<- as.factor(df.prot$Response)


data.nocb<-read_tsv("rna_assoc/rna_cpm_new.txt",col_names = TRUE,quote = "\"")
data.nocb.coll<- data.nocb[,-c(1:7)]
rna<- data.nocb.coll
rna<- as.matrix(rna)
rownames(rna)<- data.nocb$`HGNC symbol`

sample.anno<-read.csv("rna_assoc/samples_RNA.csv", header = T)
sample.anno$Tumor.response[sample.anno$Tumor.response=="Sensitive"]<- "sensitive"
sample.anno$Tumor.response[sample.anno$Tumor.response=="Refractory"]<- "refractory"
sample.anno$RNA.ID<- as.character(sample.anno$RNA.ID)

sample.anno$RNA.ID[c(7,18,19,21,27)]<- c("R69", "R80", "R82", "R85", "R93")

sample.anno.sub<- sample.anno[sample.anno$RNA.ID %in% colnames(rna),]

sample.anno.sub.ord<- sample.anno.sub[order(match(sample.anno.sub$RNA.ID, colnames(rna))),]

sen<-which(sample.anno.sub.ord$Tumor.response=="sensitive")

ref<-which(sample.anno.sub.ord$Tumor.response=="refractory")

senref<-rep(NA,106)

senref[sen]<-rep(1,length(sen)) ## 58
senref[ref]<-rep(0,length(ref)) ## 48
table(senref)

data.rna<- rna[which(rownames(rna) %in% gene),]  
df.rna<- cbind.data.frame(t(data.rna),"Response"=senref)
df.rna$Response<- as.factor(df.rna$Response)

load("CNV_coll.RData")
data.cnv<- cnv2[which(rownames(cnv2) %in% gene),]  
df.cnv2<- cbind.data.frame("cnv_id"=colnames(cnv2),t(data.cnv))
df.prot2<- cbind.data.frame("prot_id"=sample.anno.worep$Data.Label.2..data.set._.TMTplex._.sample.source._.patient.ID._sample.ID.,"Response"=senref)

dna<-read.table("purity_DNA.txt", sep="\t", header=T)
colnames(dna)[1]<- "prot_id"
colnames(dna)[2]<- "cnv_id"
dna.sub<- dna[,1:2]
dna.sub$cnv_id<- as.character(as.factor(dna.sub$cnv_id))
dna.sub$cnv_id[c(20,67,74,107)]<- c("D110","D112","D108","D89")


m<- merge(dna.sub, df.cnv2, by="cnv_id")
df.cnv<- merge(m,df.prot2,by="prot_id")
df.cnv$Response<- as.factor(df.cnv$Response)

df.cnv1<- df.cnv[,-c(1:2)]

df.prot$data<- "Protein"
df.rna$data<- "RNA"
df.cnv1$data<- "CNV"

df.prot$Response<- ifelse(df.prot$Response=="0","Refractory","Sensitive")
df.rna$Response<- ifelse(df.rna$Response=="0","Refractory","Sensitive")
df.cnv1$Response<- ifelse(df.cnv1$Response=="0","Refractory","Sensitive")

df.prot$Response<- factor(df.prot$Response, levels=c("Sensitive", "Refractory"))
df.rna$Response<- factor(df.rna$Response, levels=c("Sensitive", "Refractory"))
df.cnv1$Response<- factor(df.cnv1$Response, levels=c("Sensitive", "Refractory"))

df.cnv1$title<- "CNV" ### change to df.prot and df.rna
p.cnv<- ggplot(df.cnv1, aes(x=Response, y=TAP1, color=Response)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Response), size = 2, shape = 21, position = position_jitterdodge()) + theme_pubr() + ylab("") + scale_color_manual(values=c("#EFC000", "#0073C2")) + facet_grid(. ~ title) + theme(strip.text.x=element_text(size=30)) + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=25,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=25,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 

#### change p.cnv1 to p.rna with input df.rna and to p.prot with input df.prot

############# CPTAC-2016 Protein ##############

data<- read.table("CPTAC2/CPTAC2_coll.txt", header = T)
data.gene<- data[which(rownames(data)=="BCL2L1"),]

class<- read.table("CPTAC2_clustering.txt", header = T)

class.ord<- class[order(match(class$Sample, colnames(data))),]

identical(as.character(as.factor(class.ord$Sample)), colnames(data)) # TRUE

df.prot<- cbind.data.frame(t(data.gene),"Response"=class.ord$Response)
df.prot$Response<- as.factor(df.prot$Response)
df.prot$Response<- factor(df.prot$Response, levels=c("Sensitive", "Refractory"))

df.prot$title<- "Protein"
p.protc<- ggplot(df.prot, aes(x=Response, y=BCL2L1, color=Response)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge()) + theme_pubr() + ylab("") + scale_color_manual(values=c("#E7B800", "#0073C2")) + facet_grid(. ~ title) + theme(strip.text.x=element_text(size=30)) + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=25,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=25,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 

ggarrange(p.prot,p.rna,p.cnv, p.protc, nrow=1, ncol=4, common.legend = T, legend = "bottom")
