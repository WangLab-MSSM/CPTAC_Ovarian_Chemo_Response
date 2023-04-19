source("balanced_folds.R")
library(randomForest)
library(glmnet)
library(xgboost)
library(pROC)
library(ggplot2)
library(ggpubr)

load("meta_genes.RData")
com<- meta
load("data_prediction_FFPE_UM_primmets.RData")
load("senref_FFPE_UM_primmets.RData")

data.nocb.rib<-data.nocb.coll.nopool[which(genes %in% com),]
rownames(data.nocb.rib)<- make.names(rownames(data.nocb.rib), unique = T)

######### random Forest #########
mat.rand<-list()
for(j in 1:100)
{
  print(j)
  set.seed(j) #73, #48
  
  cv.fold=sample(rep(1:5, c(27,27,27,27,27)), 135, replace=FALSE)
  
  table(cv.fold, senref)
  
  mat<- list()
  y.pred=rep(NA, length(senref))
  for(i in 1:5)
  {
    cur.train=which(cv.fold!=i)
    cur.test=which(cv.fold==i)
    
    x.train<- data.nocb.rib[,cur.train]
    x.test<- data.nocb.rib[,cur.test]
    
    y<-as.factor(as.character(as.numeric(senref)))
    
    temp<-randomForest(y[cur.train] ~ ., data=t(x.train), importance=TRUE,proximity=TRUE, ntree=100)
    
    mat[[i]] <- temp$importance
  }
  
  mat.rand[[j]]<- mat
}

save(mat.rand,file="rf.RData")

###################### ElasticNet ######################

var.rep<-list()

for(j in 1:100)
{
  print(j)
  set.seed(j) #73, #48
  
  cv.fold=sample(rep(1:5, c(27,27,27,27,27)), 135, replace=FALSE)
  
  #cv.fold <-balanced.folds.vec(senref,5)
  
  table(cv.fold, senref)
  
  y.pred=rep(NA, length(senref))
  var<-list()
  for(i in 1:5)
  {
    cur.train=which(cv.fold!=i)
    cur.test=which(cv.fold==i)
    
    x.train=data.nocb.rib[,cur.train]
    x.test=data.nocb.rib[,cur.test]
    
    temp<-cv.glmnet(t(x.train),senref[cur.train],family = "binomial",alpha=0.1)
    var1<-coef(temp,s=temp$lambda.min)
    var[[i]]<-rownames(var1)[var1[,1]!=0][-1]
    
    y.pred[cur.test]=predict(temp, newx=t(x.test),type="response",s="lambda.min")
  }
  var.rep[[j]]<-unlist(var)
}

rib.sel<-unlist(var.rep)

table_glmnet<-cbind.data.frame("Features"=names(table(rib.sel)),"Importance"=table(rib.sel)/500)[,-2]

save(table_glmnet,file="glmnet.RData")

############# xgboost ################
library(xgboost)
library(DiagrammeR)


mat.mat<-list()
for(j in 1:100)
{
  set.seed(j) #73, #48
  
  cv.fold=sample(rep(1:5, c(27,27,27,27,27)), 135, replace=FALSE)
  #cv.fold <-balanced.folds.vec(senref,5)
  
  table(cv.fold, senref)
  
  mat<-list()
  y.pred=rep(NA, length(senref))
  for(i in 1:5)
  {
    cur.train=which(cv.fold!=i)
    cur.test=which(cv.fold==i)
    
    xtrain<-t(data.nocb.rib)[cur.train,]
    xtest<-t(data.nocb.rib)[cur.test,]
    
    temp <- xgboost(data = xtrain, label = senref[cur.train], max_depth = 2, eta = 0.1, nthread = 2, nrounds = 100, objective = "binary:logistic")
    
    y.pred[cur.test]=predict(temp, newdata=xtest)
    
    mat[[i]] <- xgb.importance(feature_names = colnames(t(data.nocb.rib)),model = temp)
  }
  mat.mat[[j]]<-mat
}

save(mat.mat,file="xgboost.RData")


fn <- function(mat.rand, th)
{
  foo <- lapply(mat.rand,function(x) do.call(rbind,x))
  foo2 <- do.call(rbind,foo)
  
  #foo2$GainFeature <- ifelse(foo2$Gain > th, foo2$Feature, NA)
  
  features<-rownames(foo2)
  rownames(foo2)<-NULL
  foo2<-data.frame(foo2)
  foo2$Feature<- features
  
  foo2$GainFeature <- ifelse(foo2$MeanDecreaseGini > th, foo2$Feature, NA)
  
  result <- as.matrix(table(foo2$GainFeature, useNA = "no"))
  result[,1] <- result[,1]/500
  
  res<- result[result[,1]>.4,]
  return(res)
}

out.xgb<-fn(mat.mat,th=.000001)
out.xgb
length(out.xgb)

out.randfor<-fn(mat.rand,th=.9)
out.randfor
length(out.randfor)

out.xgb<- as.matrix(out.xgb)
table_xgb<-cbind.data.frame(rownames(out.xgb),out.xgb[,1])
rownames(table_xgb)<- NULL
names(table_xgb)<-c("gene", "Frequency")
plot<- table_xgb[order(table_xgb$Frequency,decreasing = T),]

xg<- plot

out.randfor<- as.matrix(out.randfor)
table_randfor<-cbind.data.frame(rownames(out.randfor),out.randfor[,1])
rownames(table_randfor)<- NULL
names(table_randfor)<-c("gene", "Frequency")
plot<- table_randfor[order(table_randfor$Frequency,decreasing = T),]

rf<-plot

table_glmnet1<-table_glmnet[table_glmnet$Importance.Freq>0.4,]
names(table_glmnet1)<-c("gene", "Frequency")

plot<- table_glmnet1[order(table_glmnet1$Frequency,decreasing = T),]

glm<-plot

union(glm$gene, union(rf$gene,xg$gene))

