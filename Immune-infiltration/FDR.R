set.seed(1)
cluster<-3

list<-list.files(path = paste0("Cluster",cluster,"/"), 
           pattern = "iJRF_results_perm", all.files = TRUE,
           full.names = TRUE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = TRUE)

imp.perm.2<-imp.perm.1<-NULL
for (kk in 1:length(list)){
  load(list[kk])
  imp.perm.1<-cbind(imp.perm.1,out.new[,1])
  imp.perm.2<-cbind(imp.perm.2,out.new[,2])
  print(kk)
}

load(paste0("iJRF_results.rda"))
importance<-out.new

# -- FDR
set.seed(1)
FDR=matrix(0,50000,6)
Net=1
if (Net==1) perm=imp.perm.1
if (Net==2) perm=imp.perm.2
I<-sort(importance[,Net+2],decreasing=T)
TH=1E-4

for (j in 10:30000){
  FDR[j,Net]<-((sum(perm>I[j])/dim(perm)[2]))/j
  if (FDR[j,Net]>TH) break
  print(c(j,FDR[j,Net]))
}

index<-grep("Location",importance[,1])
importance<-importance[-index,]

index<-grep("Location",importance[,2])
importance<-importance[-index,]

importance<-importance[importance[,Net+2]>I[j],seq(1,2)]


if (Net==1) type="R"
if (Net==2) type="S"
file=paste0("Network_global_type_",type,"_TH_",TH,".txt")

write.table(importance,file=file,
            row.names=FALSE,quote=FALSE,append=FALSE,sep="\t",col.names=c("gene1","gene2"))
