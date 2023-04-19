fin64<- read.csv("Final_64_proteins.csv")
com<- fin64$Proteins

######### FFPE Discovery and FFPE Pilot ##########
load("data_prediction_FFPE_UM_primmets_noneo.RData") #### Discovery
load("senref_FFPE_UM_primmets_noneo.RData")

load("data_prediction_UM_Pilot.RData") #### Pilot
load("senref_UM_Pilot.RData")

data.nocb.rib<-data.nocb.coll.nopool[which(rownames(data.nocb.coll.nopool) %in% com),]

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

roc.all1 <- roc(senref.pool1, y.pred.avg.pool1.all, ci = TRUE, direction="<")

ciobj1<-ci.se(roc.all1, specificities = seq(0,1,.01))
dat.ci1 <- data.frame(x = as.numeric(rownames(ciobj1)),
                      lower = ciobj1[, 1],
                      upper = ciobj1[, 3])

######### Frozen ##########

load("data_nov.RData"); load("senref_nov.RData")
load("data_prediction_FFPE_UM_primmets_noneo.RData")
load("senref_FFPE_UM_primmets_noneo.RData")

######### For both non-overlapping and overlapping samples in Frozen 
data.nocb.rib<-data.nocb.coll.nopool1[which(rownames(data.nocb.coll.nopool1) %in% com),]
rownames(data.nocb.rib)<- make.names(rownames(data.nocb.rib), unique = T)

data.nocb.rib.FFPE<- data.nocb.rib

data.nocb.rib.FFPE1<- data.nocb.rib.FFPE[order(match(rownames(data.nocb.rib.FFPE),rownames(data.nocb.rib.Frozen))),]

senref.FZ<- senref
senref.FFPE<- senref1

x.train<- data.nocb.rib.FFPE1
x.test<- data.nocb.rib.Frozen

set.seed(1000000)
temp.en<-try(cv.glmnet(t(x.train),senref.FFPE,family = "binomial", alpha=0.2, nfolds = 3), silent = T)
coef<-coef(temp.en, s=c(temp.en$lambda.min,temp.en$lambda.1se))

y.pred.EN=predict(temp.en, newx=t(x.test),type="response",s="lambda.min")
roc.EN <- roc(senref.FZ, as.numeric(y.pred.EN), ci = TRUE, direction="<")

y.pred.EN.FFPE=predict(temp.en, newx=t(x.train),type="response",s="lambda.min")

temp.xg <- xgboost(data = t(x.train), label = senref.FFPE, max_depth = 2, eta = 0.1, nthread = 2, nrounds = 500, objective = "binary:logistic")

y.pred.XG=predict(temp.xg, newdata=t(x.test), ntreelimit = 500)
roc.XG <- roc(senref.FZ, y.pred.XG, ci = TRUE, direction="<")

y.pred.XG.FFPE=predict(temp.xg, newdata=t(x.train), ntreelimit = 500)

y<-as.factor(as.character(as.numeric(senref.FFPE)))

temp.rf<-randomForest(y ~ ., data=t(x.train), importance=TRUE,proximity=TRUE, ntree=1000)
y.pred.RF=predict(temp.rf, newdata=t(x.test),type="prob")[,2]

y.pred.RF.FFPE=predict(temp.rf, newdata=t(x.train),type="prob")[,2]

y.pred.avg1.all<- apply(cbind(y.pred.XG, y.pred.RF, y.pred.EN),1,mean)

roc.avg.all <- roc(senref.FZ, y.pred.avg1.all, ci = TRUE, direction="<")


###### roc plot #########
gtrain<-roc.all1
ggroc(gtrain,alpha = 0.5, colour = "green4", linetype = 1, size = 2) +  theme_bw() + annotate(geom="text", x=0.2, y=0.4, label="AUC = 0.83", color="red",size=7) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))  + geom_ribbon(data = dat.ci1, aes(x = x, ymin = lower, ymax = upper), fill = "green", alpha= 0.2) 