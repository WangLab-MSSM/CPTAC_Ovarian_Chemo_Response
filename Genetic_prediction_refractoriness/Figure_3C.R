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
location[ov]<-rep("ov",length(ov)) ## 71
location[om]<-rep("om",length(om)) ## 51

df.prot1<- cbind.data.frame("prot_id"=sample.anno.worep$Data.Label.2..data.set._.TMTplex._.sample.source._.patient.ID._sample.ID., "loc"=location, "age"=age)

loh<- read.csv("chr17LOHdata.csv", header=T)

loh1<- loh[-c(8,11,13),]
colnames(loh1)[1]<- "cnv_id"

mut<- read.csv("brca_relaxed_mutations.csv", header = T)
mut1<- mut[-c(7,10,121),]
colnames(mut1)[2]<- "cnv_id"

dna<-read.table("purity_DNA.txt", sep="\t", header=T)
colnames(dna)[1]<- "prot_id"
colnames(dna)[2]<- "cnv_id"
dna.sub<- dna[,1:2]
dna.sub$cnv_id<- as.character(as.factor(dna.sub$cnv_id))
dna.sub$cnv_id[c(20,67,74,107)]<- c("D110","D112","D108","D89")

m<- merge(dna.sub,loh1,by="cnv_id")
mm<- merge(m, mut1, by="cnv_id")

load("data_prediction_FFPE_UM_primmets_noneo.RData")
load("senref_FFPE_UM_primmets_noneo.RData")

df.prot<- cbind.data.frame("prot_id"=colnames(data.nocb.coll.nopool1), "senref"=senref1)

df.prot2<- merge(df.prot,df.prot1, by="prot_id")

mmm<- merge(mm, df.prot2, by="prot_id")

data.nocb.coll.nopool2<- data.nocb.coll.nopool1[,which(colnames(data.nocb.coll.nopool1) %in% mmm$prot_id)]

mmm.ord<- mmm[order(match(mmm$prot_id, colnames(data.nocb.coll.nopool2))),]

senref<- mmm.ord$senref

brca<- ifelse(mmm.ord$brca_mut == "TRUE", 1,0)

###### for clin var + brca mut + chr17LOH ######
data.nocb.rib<- rbind.data.frame("brca"=brca,"loc"=mmm.ord$loc,"age"=mmm.ord$age, "loh"=mmm.ord$chr17LOH)

###### for clin var + brca mut ######
data.nocb.rib<- rbind.data.frame("brca"=brca,"loc"=mmm.ord$loc,"age"=mmm.ord$age)


source("balanced_folds.R")

y.pred.pool.XG <- NULL
y.pred.pool.RF <- NULL
y.pred.pool.EN <- NULL

y.pred.avg.pool.all <- NULL

y.pred.avg.pool1.all <- NULL

for(j in 1:100)
{
  print(j)
  set.seed(j)
  cv.fold <-balanced.folds.vec(senref,5)
  
  table(cv.fold, senref)
  
  y.pred.EN=rep(NA, length(senref))
  y.pred.XG=rep(NA, length(senref))
  y.pred.RF=rep(NA, length(senref))
  
  y.pred.avg.all<- rep(NA, length(senref))
  
  y.pred.XG1=rep(NA, length(senref))
  y.pred.EN1=rep(NA, length(senref))
  y.pred.RF1=rep(NA, length(senref))
  
  
  for(i in 1:5)
  {
    cur.train=which(cv.fold!=i)
    cur.test=which(cv.fold==i)
    
    ##### EN #######
    x.train.en<- data.nocb.rib[,cur.train]
    x.test.en<- data.nocb.rib[,cur.test]
    
    temp.en<-try(cv.glmnet(t(x.train.en),senref[cur.train],family = "binomial", alpha=0.2, nfolds = 3), silent = T)
    
    if(class(temp.en)=="try-error"){
      stop = TRUE # Fire the flag, and break the inner loop
      break;
    }
    y.pred.EN[cur.test]=predict(temp.en, newx=t(x.test.en),type="response",s="lambda.min")
    
    y.pred.EN1[cur.test]=predict(temp.en, newx=t(x.test.en),type="response",s="lambda.min")
    y.pred.EN1[cur.train]=predict(temp.en, newx = t(x.train.en),type="response",s="lambda.min")
    
    ######################## xgb ############################
    x.train.xg<- data.nocb.rib[,cur.train]
    x.test.xg<- data.nocb.rib[,cur.test]
    
    temp.xg <- xgboost(data = t(x.train.xg), label = senref[cur.train], max_depth = 2, eta = 0.1, nthread = 2, nrounds = 100, objective = "binary:logistic")
    
    y.pred.XG[cur.test]=predict(temp.xg, newdata=t(x.test.xg))
    
    y.pred.XG1[cur.test]=predict(temp.xg, newdata=t(x.test.xg))
    y.pred.XG1[cur.train]=predict(temp.xg, newdata=t(x.train.xg))
    
    ################################ RF ###############################
    x.train.rf<- data.nocb.rib[,cur.train]
    x.test.rf<- data.nocb.rib[,cur.test]
    
    y<-as.factor(as.character(as.numeric(senref)))
    
    set.seed(1000000)
    temp.rf<-randomForest(y[cur.train] ~ ., data=t(x.train.rf), importance=TRUE,proximity=TRUE, ntree=100)
    y.pred.RF[cur.test]=predict(temp.rf, newdata=t(x.test.rf),type="prob")[,2]
    
    
    y.pred.RF1[cur.test]=predict(temp.rf, newdata=t(x.test.rf),type="prob")[,2]
    y.pred.RF1[cur.train]=predict(temp.rf, newdata=t(x.train.rf),type="prob")[,2]
    
    ############################ aggregation ###############################
    ########################################################################
    ##################### XG ###############################
    
    ############################## all ###########################    
    
    y.pred.train.all <- cbind.data.frame(senref = senref[cur.train], en = y.pred.EN1[cur.train], xg = y.pred.XG1[cur.train], rf = y.pred.RF1[cur.train])
    
    
    y.pred.mat.all<-cbind.data.frame(en = y.pred.EN1[cur.test], xg = y.pred.XG1[cur.test], rf = y.pred.RF1[cur.test])
    
    temp.avg.all<- glm(senref ~  en + xg + rf, family=binomial(link = "logit"), data = y.pred.train.all)
    
    y.pred.avg.all[cur.test]<- predict(temp.avg.all, newdata=y.pred.mat.all, type="response")
  }  
  
  ######################## predict: simple avg of 3 pathways ######################
  
  y.pred.mat1.all<-cbind(y.pred.EN, y.pred.XG,y.pred.RF)
  
  y.pred.avg1.all<-apply(y.pred.mat1.all,1,mean)
  
  ####################### pooling indiv methods #####################
  y.pred.pool.XG<- c(y.pred.pool.XG,y.pred.XG)
  y.pred.pool.RF<- c(y.pred.pool.RF,y.pred.RF)
  y.pred.pool.EN<- c(y.pred.pool.EN,y.pred.EN)
  
  ############# pooling weighted avg of pathways ####################    
  
  y.pred.avg.pool.all<- c(y.pred.avg.pool.all,y.pred.avg.all)
  
  ############# pooling simple avg of pathways ####################    
  
  y.pred.avg.pool1.all<- c(y.pred.avg.pool1.all,y.pred.avg1.all)
}

y.pred.avg.pool1.all<-apply(matrix(y.pred.avg.pool1.all,length(senref),100),1,function(x) mean(x, na.rm=T))

senref.pool1<- senref

roc.clin.brca.chr17loh <- roc(senref.pool1, y.pred.avg.pool1.all, ci = TRUE, direction="<")
roc.clin.brca <- roc(senref.pool1, y.pred.avg.pool1.all, ci = TRUE, direction="<")

roc.list<- list("BRCA + clin var"=roc.clin.brca, "BRCA + clin var + chr17LOH"=roc.clin.brca.chr17loh)

ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, 0.01))

dat.ci.list <- lapply(ci.list, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

ggroc(roc.list, aes("color", "linetype"), size = 2) +
  annotate(geom="text", x=0.3, y=0.25, label="AUC = 0.61", color="#F8766D",size=7) +
  
  annotate(geom="text", x=0.3, y=0.2, label="AUC = 0.66", color="#009900",size=7) +
  
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid") +
  
  geom_ribbon(data = dat.ci.list[[1]], aes(x = x, ymin = lower, ymax = upper), fill = "red", alpha= 0.2,
              inherit.aes = F) +
  
  geom_ribbon(data = dat.ci.list[[2]], aes(x = x, ymin = lower, ymax = upper), fill = "green", alpha= 0.2,
              inherit.aes = F) + theme_bw() +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain")) 
