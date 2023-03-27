
  # ---- Input   -------------------------------------------------------------------------------------------- ###

  #install.packages("iJRFNet_1.1-4.tar.gz")
  library(iJRFNet)
  library(stringr)
  
  cluster=1
  perm.seed=1
  
  ### --- Load Data ----------------------------------------------------------------------------------------- ###
  load(paste0("Data_cluster_",cluster,".rda"))
  geneID<-rownames(data.R)
  
  # -- Run iJRFNet
  dummy<-function(location,data){
    data<-t(apply(data,1,as.numeric))
    data<-t(apply(data,1,function(x)(x-mean(x))/sd(x)))
    
    location.dummy<-matrix(0,length(location),3)
    location.dummy[location=="OM",1]<-1
    location.dummy[location=="OV",2]<-1
    location.dummy[location=="Other",3]<-1
    colnames(location.dummy)<-c("Location-OM","Location-OV","Location-Other")
    data<-rbind(t(location.dummy),data)
    return(data)
  }
  data.R<-dummy(location.R,data.R)
  data.S<-dummy(location.S,data.S)
  
  out.new<-iJRFNet_permutation(X=list(data.R,data.S),genes.name=c(rownames(data.R)), model="iJRF",seed=perm.seed,to.store=100000)

    save(iJRF_results_perm_",perm.seed,".rda",sep="")) # -- save on file importance scores
