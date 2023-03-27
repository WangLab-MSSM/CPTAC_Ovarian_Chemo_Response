library("preprocessCore")
source("repbulk_unconstrained_weights_slab_missing_samples.R")
library(stringr)

# -- load global proteomic
data<-read.table("FD_GLBLprot_MI_FDbridge_Abund_20201002_Imput_v01.tsv",sep="\t",header=TRUE)
geneID<-data[,2]
data<-data[,-seq(1,9)]
sample<-colnames(data)
data<-apply(data,2,as.numeric)
rownames(data)<-geneID

# -- load RNA
RNA<-read.table("RNA.tsv",header=TRUE,row.names=1)
sample.rna<-colnames(RNA)

# -- map IDs
clinical<-read.table("Clinical_TMT_order.csv",header=TRUE,sep=",",stringsAsFactors =FALSE)
mg<-match(str_replace(str_replace(str_replace(str_replace(str_replace(sample,"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),str_replace(str_replace(str_replace(str_replace(str_replace(clinical[,1],"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"))
clinical<-clinical[mg,]

mg<-match(clinical[,"RNA.ID"],sample.rna)
RNA.final<-matrix(NA,dim(RNA)[1],dim(clinical)[1])
RNA.final[,!is.na(mg)]<-as.matrix(RNA[,mg[!is.na(mg)]])
rownames(RNA.final)<-rownames(RNA)

# -- load scRNA signature
sc.immune<-read.table("~/Box Sync/PTRC/Signature_matrix/SC_GeneList.csv",sep=",",fill=TRUE,header=TRUE,stringsAsFactors =FALSE)

# -- load LM22
LM22<-read.table("LM22.csv",sep=",",fill=TRUE,header=TRUE,row.names=1)
LM22<-cbind(LM22[,-seq(17,18)],rowSums(LM22[,seq(17,18)])/2) # dendritic
LM22<-cbind(LM22[,-seq(11,12)],rowSums(LM22[,seq(11,12)])/2) # nk cells
LM22<-cbind(LM22[,-seq(15,16)],rowSums(LM22[,seq(15,16)])/2) # mast cell
LM22<-cbind(LM22[,-seq(1,2)],rowSums(LM22[,seq(1,2)])/2) # b cell

colnames(LM22)<-c(colnames(LM22)[c(seq(1,14))],"Dendritic","NK Cells","Mast Cells","B Cells")
cell.type<-c(colnames(LM22),colnames(sc.immune))
lm22<-dim(LM22)[2]

gene.single<-c(as.matrix(sc.immune)); gene.single<-gene.single[gene.single!=""]
mg<-match(rownames(data),unique(c(rownames(LM22),gene.single)))
data<-data[!is.na(mg),]

mg<-match(rownames(RNA.final),unique(c(rownames(LM22),gene.single)))
RNA.final<-RNA.final[!is.na(mg),]

# --- build index matrix
build.matrix<-function(data,LM22,cell.type,lm22){
index.matrix<-c(NULL,NULL,NULL)
for (s in 1:length(cell.type)){
  for (k in 1:length(cell.type)) {
    
    if (k!=s){
      if (s <= lm22 & k <= lm22){
        i<-sort(LM22[,s],decreasing=TRUE,index=TRUE)
        
        gene<-NULL
        sss<-0
        for (ss in 1:dim(LM22)[1]){
          compare<-LM22[i$ix[ss],c(s,k)] 
          if (compare[1]>5*compare[2]) {
            gene<-c(gene,rownames(LM22)[i$ix[ss]])
            sss<-sss+1
            if (sss==10) break;
            if (ss==dim(LM22)[1]) print("yes")
          }
        }
        gene.s<-gene
        index.s<-cbind(rep(s,length(gene.s)),rep(k,length(gene.s)),gene.s)
        
      }
      
      if (s <= lm22 & k > lm22){
        
        i<-sort(LM22[,s],decreasing=TRUE,index=TRUE)
        gene<-rownames(LM22)[i$ix[seq(1,30)]]
        gene<-gene[is.na(match(gene,sc.immune[,k-lm22]))]
        gene<-gene[seq(1,10)]
        index.s<-cbind(rep(s,length(gene)),rep(k,length(gene)),gene)
        
      }
      
      if (k <= lm22 & s > lm22){
        #      i<-sort(LM22[,k],decreasing=TRUE,index=TRUE)
        #      gene<-rownames(LM22)[i$ix[seq(1,20)]]
        gene<-sc.immune[,s-lm22]
        gene<-gene[!is.na(match(gene,rownames(data)))]
        index.s<-cbind(rep(s,length(gene)),rep(k,length(gene)),gene)
      }
      
      if (k > lm22 & s > lm22){
        gene<-sc.immune[,s-lm22]
        gene<-gene[is.na(match(gene,sc.immune[,k-lm22]))]
        gene<-gene[!is.na(match(gene,rownames(data)))]
        if (length(gene)>0) index.s<-cbind(rep(s,length(gene)),rep(k,length(gene)),gene)
      }
      index.matrix<-rbind(index.matrix,index.s)
    }
  }
}
index.matrix<-index.matrix[index.matrix[,3]!="",]
return(index.matrix)
}

LM22.protein<-LM22[!is.na(match(rownames(LM22),rownames(data))),]
index.matrix.pro<-build.matrix(data,LM22.protein,cell.type,lm22)

LM22.rna<-LM22[!is.na(match(rownames(LM22),rownames(RNA.final))),]
index.matrix.rna<-build.matrix(RNA.final,LM22.rna,cell.type,lm22)

mg<-match(index.matrix.pro[,3],rownames(data))
index.matrix.pro[,3]<-mg[!is.na(mg)]

mg<-match(index.matrix.rna[,3],rownames(RNA.final))
index.matrix.rna[,3]<-mg[!is.na(mg)]+dim(data)[1]

# -- final data
data<-rbind(data,RNA.final)
index.matrix<-rbind(index.matrix.pro,index.matrix.rna)
index.matrix<-apply(index.matrix,2,as.numeric)

n.iter<-2000
burn.in<-1000
k.fix<-length(cell.type)
data<-t(apply(data,1,function(x)(x-mean(x[!is.na(x)]))/sd(x[!is.na(x)])))

gibbs<-gibbs.sampling(n.iter=n.iter,data,p=dim(data)[1],n=dim(data)[2],k.fix,index.matrix,burn.in,mean.prior=matrix(0,dim(data)[1],k.fix),sigma.prior=1)

pi.post<-matrix(0,dim(data)[2],k.fix)
mu.post<-matrix(0,dim(data)[1],k.fix)
for (k in 1:k.fix) {
  pi.post[,k]<-apply(gibbs[[1]][[k]],1,median)
  mu.post[,k]<-apply(gibbs[[2]][[k]],1,median)
}


library(ggpubr)
colnames(pi.post)<-colnames(mu.post)<-c(cell.type)
rownames(pi.post)<-colnames(data)
rownames(mu.post)<-rownames(data)

save(gibbs,pi.post,mu.post,cell.type,file=paste0("results.rda"))

