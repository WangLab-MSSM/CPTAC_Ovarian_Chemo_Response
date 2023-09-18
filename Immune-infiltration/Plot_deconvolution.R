library(gginnards)
library(stringr)
library(ggpubr)
library("ComplexHeatmap")
library(RColorBrewer) 
library("circlize")

source("function/functions.R")

# -- load data
data<-read.table("RNA.tsv",header=TRUE,row.names=1)
geneID<-rownames(data)
sample<-colnames(data)
data<-apply(data,2,as.numeric)
rownames(data)<-geneID

# -- load pathway scores
hallmark<-read.table("data/ssGSEA/Hallmark.tsv",header=TRUE,fill=TRUE)

# --- load BayesDeBulk deconvolution results
load("results/results.rda")
for (s in 1:dim(pi.post)[1]) pi.post[s,]<- pi.post[s,]/sum( pi.post[s,]) 

mg<-match(rownames(pi.post),colnames(hallmark))
pi.post<-pi.post[!is.na(mg),]; hallmark<-hallmark[,mg[!is.na(mg)]]

# -- load factor model deconvolution
cell.factor<-read.table("results/Factor_model_deconvolution.txt",header=TRUE,row.names=1,sep="\t")
mg<-match(map(rownames(cell.factor)),map(rownames(pi.post)))
cell.factor<-cell.factor[!is.na(mg),]
pi.post<-pi.post[mg[!is.na(mg)],]
hallmark<-hallmark[,mg[!is.na(mg)]]
mean.post.global<-pi.post

# -- lad protein clusters
cluster<-read.table("results/FFPE_protein_subtypes.csv",header=FALSE,row.names=1,sep=",")
mg<-match(map(rownames(cell.factor)),map(rownames(cluster)))
cluster<-cluster[mg,1]

# -- load protein data
global<-read.table("data/FD_GLBLprot_MI_FDbridge_Abund_20201002_Imput_v01.tsv",header=TRUE,row.names=1,sep="\t")
rownames(global)<-global[,4]
mg<-match(map(rownames(cell.factor)),map(colnames(global)))
global<-global[,mg]

global<-global[c("CD8A","CD4","STAT1","CXCL10","IDO1", "HLA-A" ,   "HLA-B"  ,  "HLA-C","HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2" ),]

# -- Load pd1 signature
pd1<-read.table(paste0("results/PD1_proteomic_score.csv"),sep=",",header=TRUE)
mg<-match(str_replace(str_replace(str_replace(str_replace(str_replace(rownames(mean.post.global),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),str_replace(str_replace(str_replace(str_replace(str_replace(colnames(pd1),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"))
pd1<-pd1[1,mg]

# -- load clinical data
library(stringr)
clinical<-read.table("data/Clinical_TMT_order.csv",header=TRUE,sep=",",stringsAsFactors =FALSE)

clinical.sample<-read.table("data/sampleID_annotation.csv",header=TRUE,sep=",",stringsAsFactors =FALSE)
mg<-match(str_replace(str_replace(str_replace(str_replace(str_replace(rownames(mean.post.global),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),str_replace(str_replace(str_replace(str_replace(str_replace(clinical[,1],"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"),"[[.]]","-"))
clinical<-clinical[mg,]

mean.post.global<-mean.post.global[,apply(mean.post.global,2,sd)!=0]
mean.post.global<-mean.post.global[,colSums(mean.post.global>0.05)>.3*dim(mean.post.global)[1]]
mean.post.global<-t(apply(mean.post.global,2,function(x)(x-mean(x))/sd(x)))
cell.factor<-t(apply(cell.factor,2,function(x)(x-mean(x))/sd(x)))

response<-as.character(clinical[,"Tumor.response"])
response[response=="sensitive"]<-"Sensitive"; response[response=="refractory"]<-"Refractory"

location<-as.character(clinical[,"Tumor.Location.Group"])
location[location!="OV" & location!="OM"]<-"Other"
ha<-HeatmapAnnotation(df=data.frame(Anti.PD1=as.numeric(pd1),
  Stage=clinical[,"Tumor.Stage"],  
                                    Tumor.location=as.character(location),
  tumor.response=paste(as.character(clinical[,"Tumor.type"]),as.character(clinical[,"Tumor.response"])),
  Cluster=as.character(cluster),
  tumor=as.character(clinical[,"Tumor.type"]),
  response=as.character(clinical[,"Tumor.response"])),
  col = list( Cluster=c( "2"= "#B30303", "3"="gold",   "4"="darkblue",   "5"="lightblue","1"="darkgreen"),
                Stage=c("I"="gray95","II"="gray85","III"="gray65","IV"="gray30"),
                Tumor.location=c( "OV"= "#B30303", "OM"="darkblue",   "Other"="gray80"),
              tumor.response=c( "Metastatic refractory"= "#B30303", "Metastatic sensitive"="gold",   "Primary refractory"="darkblue",   "Primary sensitive"="lightblue"),
              response=c("sensitive" = "gray80", 
                         "refractory" = "black"),
              tumor=c("Metastatic"="black","Primary"="gray80","Mix"="gold")),
  na_col="gray97",
  show_annotation_name = TRUE,annotation_name_gp=gpar(fontsize = 9,fontface="bold"),
  simple_anno_size = unit(.35, "cm"),gap=unit(.11, "cm"))

mean.post.global<-mean.post.global[rownames(mean.post.global)!="Neutrophils",]
m1<-hclust(dist(mean.post.global))$order

hallmark<-t(apply(hallmark,1,function(x)(x-mean(x))/sd(x)))
pathway<-c("HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_ANGIOGENESIS","HALLMARK_COMPLEMENT","HALLMARK_G2M_CHECKPOINT" ,"HALLMARK_DNA_REPAIR" ,
           "HALLMARK_INFLAMMATORY_RESPONSE" ,"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_OXIDATIVE_PHOSPHORYLATION" ,"HALLMARK_IL2_STAT5_SIGNALING" ,"HALLMARK_TGF_BETA_SIGNALING" ,"HALLMARK_MYC_TARGETS_V1","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_GLYCOLYSIS","HALLMARK_E2F_TARGETS","HALLMARK_FATTY_ACID_METABOLISM"  )

pathway.plot<-as.matrix(hallmark[pathway,])
rownames(pathway.plot)<- str_to_title(str_replace(str_replace(str_replace(str_replace(str_replace(rownames(pathway.plot),"HALLMARK_",""),"_"," "),"_"," "),"_"," "),"_"," s"))

global<-t(apply(global,1,function(x)(x-mean(x))/sd(x)))
row.split<-data.frame(factor(as.character(c(1,2,2,2,1,2,2,2,2,2,2,4,3,1,4,2,1,4,3,5,5,5,5,5,5,5,5,5,5,5,5,rep(6,length(pathway)))[-c(1,3,4,5)]),levels=c("1","3","4","6","2","5")))
ht1=Heatmap(rbind(as.matrix(mean.post.global[-c(1,3,4,5),]),as.matrix((cell.factor)),as.matrix(global),as.matrix(pathway.plot)),
             col = colorRamp2(c(seq(-3,-.8,length=100),seq(-.79,.79,length=100),seq(.8,3,length=100)), 
                               colorRampPalette(brewer.pal(9, "RdBu"))(300)[seq(300,1)]),#row_km = 4,
            row_split = row.split,
column_split=data.frame(factor(as.character(cluster)),clinical[,"Tumor.response"]),
            show_row_names=TRUE,show_column_names=FALSE,column_dend_reorder=FALSE,show_column_dend=FALSE,row_dend_reorder=FALSE,
            column_dend_height=unit(1, "mm"),show_row_dend=FALSE,
            row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 6),cluster_columns =TRUE,cluster_rows =TRUE,
            top_annotation=ha,
            heatmap_legend_param = list(title = " ",legend_direction = "vertical",
                                        legend_width = unit(6, "cm"),
                                        legend_height = unit(4, "cm")))


pdf("Heatmap_deconvolution.pdf",width=11,height=7)
draw(ht1, heatmap_legend_side = "left", annotation_legend_side="left",newpage=FALSE)
dev.off()


# -- Generate boxplot
width=15
height=8
data.plot<-(mean.post.global)
response=as.character(clinical[,"Tumor.response"])

df<-data.frame(value=c(t(data.plot)),
               cell=rep(rownames(data.plot),each=dim(data.plot)[2]),
               response=rep(as.character(response),dim(data.plot)[1]),
               tumor=rep(as.character(clinical[,"Tumor.type"]),dim(data.plot)[1]),
               type=paste(rep(as.character(response),dim(data.plot)[1]),rep(as.character(clinical[,"Tumor.type"]),dim(data.plot)[1])),
               treatment=rep(as.character(clinical[,"Neo.adjuvant"]),dim(data.plot)[1]),
               location=rep(as.character(clinical[,"Tumor.Location.Group"]),dim(data.plot)[1]),
               cluster=paste("cluster",as.factor(rep(cluster,dim(data.plot)[1]))))



theme_set(theme_grey(base_size = 20))

pdf(paste0("Boxplot_protein_cluster_deconvolution.pdf"),width=width,height=height)
q<-ggplot(df,aes(cluster,value,color=response))+geom_boxplot() + geom_point(position = position_jitterdodge()) +facet_wrap(~cell,scales="free_x",ncol=6)+theme_bw() +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 0.5, vjust = -1) +
  scale_color_manual(values=c("#0073C2","#E7B800"))+
  xlab(" ")+ylab(" ")+
  theme_bw()+
  theme(text=element_text(size=9,color = "black",face = "bold"))+coord_flip()+
  theme(axis.text.x = element_text( hjust = 1,size=11,color = "black",face = "bold"))+
  theme(axis.text.y = element_text( hjust = 1,size=11,color = "black",face = "bold"))+
  theme(plot.title=element_text(size=9),strip.text = element_text(size = 12),
        axis.title = element_text( size = 11, face = "bold" ))+
  ylab(" ")+xlab(" ")+ggtitle(" ")+
  theme(plot.title = element_text(size=15))+ theme(
    axis.ticks.y=element_blank())+theme(legend.position="none")+theme(strip.background = element_blank())

print(q)
dev.off()


# -- Plot Pd1 score
df<-data.frame(value=as.numeric(pd1),cluster=paste0("Cluster ",cluster),response=response)

pdf(paste0("Boxplot_protein_cluster_PD1.pdf"),width=3.3,height=2.8)
q<-ggplot(df,aes(cluster,value,color=response))+theme_bw() +
  geom_boxplot() + geom_point(position = position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 0.5, vjust = -1) +
  scale_color_manual(values=c("#0073C2","#E7B800"))+
  xlab(" ")+ylab(" ")+
  theme_bw()+
  theme(text=element_text(size=9,color = "black",face = "bold"))+coord_flip()+
  theme(axis.text.x = element_text( hjust = 1,size=11,color = "black",face = "bold"))+
  theme(axis.text.y = element_text( hjust = 1,size=11,color = "black",face = "bold"))+
  theme(plot.title=element_text(size=9),strip.text = element_text(size = 12),
        axis.title = element_text( size = 11, face = "bold" ))+
  ylab(" ")+xlab(" ")+ggtitle(" ")+
  theme(plot.title = element_text(size=15))+ theme(
    axis.ticks.y=element_blank())+theme(legend.position="none")+theme(strip.background = element_blank())

print(q)
dev.off()
