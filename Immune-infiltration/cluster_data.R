library(stringr)

### --- Load Data ----------------------------------------------------------------------------------------- ###

data<-read.table("FD_GLBLprot_MI_FDbridge_Abund_20201002_Imput_v02.tsv",header=TRUE,row.names=1)
geneID<-as.character(data[,2])
data<-data[,-seq(1,9)]

sample<-colnames(data)

# -- load clinical data
clinical<-read.table("Clinical_TMT_order.csv",header=TRUE,sep=",",stringsAsFactors =FALSE)
mg<-match(str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(sample,"[[.]][[.]]","."),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),str_replace(str_replace(str_replace(str_replace(str_replace(clinical[,1],"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"))
clinical<-clinical[mg,]

patient.id<-clinical[,"Patient.ID"]

# -- merge multiple patient id
duplicate<-unique(patient.id[duplicated(patient.id)])

index<-(is.na(match(patient.id,duplicate)))
data.new<-data[,index]
clinical.new<-clinical[index,]
sample.new<-sample[index]
patient.new<-c(patient.id[index],duplicate)

for (s in 1:length(duplicate)){
  index<-(patient.id==duplicate[s])
  print(sum(index))
  new<-apply(data[,index],1,mean)
  data.new<-cbind(data.new,new)
}
clinical.new<-rbind(clinical.new,clinical[match(duplicate,patient.id),])

data<-data.new
clinical<-clinical.new
sampleID<-patient.new

# - k means clustering
K<-3
set.seed(1)
cluster<-kmeans((data),K)[[1]]
table(cluster)

index<-(clinical[,"Tumor.response"]=="refractory") 
data2<-data[,index]
clinical.2<-clinical[index,]

index<-(clinical[,"Tumor.response"]=="sensitive") 
data1<-data[,index]
clinical.1<-clinical[index,]

data1<-t(apply(data1,1,as.numeric))
data1<-t(apply(data1,1,function(x)(x-mean(x))/sd(x)))
rownames(data1)<-geneID

data2<-t(apply(data2,1,as.numeric))
data2<-t(apply(data2,1,function(x)(x-mean(x))/sd(x)))
rownames(data2)<-geneID

#-- save cluster data
clinical.R<-clinical.2
clinical.S<-clinical.1

location.R<-clinical.R[,"Tumor.Location.Group"]
location.R[location.R!="OV" & location.R!="OM"]<-"Other"

location.S<-clinical.S[,"Tumor.Location.Group"]
location.S[location.S!="OV" & location.S!="OM"]<-"Other"

for (s in 1:K){
  data.R<-data2[cluster==s,]
  data.S<-data1[cluster==s,]

  gene<-geneID[cluster==s]
  print(length(gene))
  save(gene,data.S,data.R,location.R,location.S,file=paste0("~/Box Sync/PTRC/Network_analysis/global/Cluster/Data_cluster_",s,".rda"))
}
