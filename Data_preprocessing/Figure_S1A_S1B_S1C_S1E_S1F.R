############ S1A #############

library(ggpubr)
library(ggplot2)

sample.anno<- read.csv("samples_FFPE_Aug19.csv")
br<- which(grepl("bridge",sample.anno$Process.ID)==T)
sample.anno<- sample.anno[-br,]

id<-c("D8","E5","H3","M8","R8")

rep1<-which(sample.anno$Process.ID=="D8") 
rep2<-which(sample.anno$Process.ID=="E5") 
rep3<-which(sample.anno$Process.ID=="H3") 
rep4<-which(sample.anno$Process.ID=="M8") 
rep5<-which(sample.anno$Process.ID=="R8") 

sample.anno.worep<-sample.anno[-c(rep1[-1],rep2[-1],rep3[-1],rep4[-1],rep5[-1]),]

cov<- read.csv("covariates.csv", header = T)
cov.ord<- cov[order(match(cov$Process.ID,sample.anno.worep$Process.ID)),]

fh<- which(sample.anno.worep$Sample.Source=="FH")
mc<- which(sample.anno.worep$Sample.Source=="MC")

site<- rep("UA",158)

site[fh]<- rep("FH",length(fh))
site[mc]<- rep("MC",length(mc))

age<- cov.ord$Sample.age..years.

df<- cbind.data.frame("site" = site, "age" = age)

ggplot(df, aes(x=site, y=age)) + geom_boxplot() + geom_jitter() + theme_pubr() + xlab("sample source site") + ylab("sample age (years)") + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.position = "none") 


############# S1B #############
data<- read.csv("FFPE_protein_batch.tsv", sep="\t")

data.nocb<-data[,-c(1:9)]


panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y,use="complete.obs"), digits=3)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}

rep1<-which(sample.anno$Process.ID=="D8") 
rep2<-which(sample.anno$Process.ID=="E5") 
rep3<-which(sample.anno$Process.ID=="H3") 
rep4<-which(sample.anno$Process.ID=="M8") 
rep5<-which(sample.anno$Process.ID=="R8") 

# Create the plots
pairs(data.nocb[,rep1], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

################ S1C ###############
df<- read.csv("FFPEFrozen_overlap_corr_protein.csv") ##### load from Table S1

hist(df$Correlation, xlab="Correlation (R)", main = "FFPE-Frozen sample pairs")

########## S1E and S1F ############
###### protein
data<- read.csv("FFPE_protein_batch.tsv", sep="\t")

data.nocb<-data[,-c(1:9)]

row.miss<- apply(data.nocb,1, function(x) mean(is.na(x)))

hist(row.miss, xlab="Missing rate per protein", main = "Distribution of missing
rate per protein")

###### phospho
data<- read.csv("FFPE_phos_batch.tsv", sep="\t")

data.nocb<-data[,-c(1:8)]

row.miss<- apply(data.nocb,1, function(x) mean(is.na(x)))

hist(row.miss, xlab="Missing rate per phosphosite", main = "Distribution of missing
rate per phosphosite")

########## S1G was generated using excel ############
