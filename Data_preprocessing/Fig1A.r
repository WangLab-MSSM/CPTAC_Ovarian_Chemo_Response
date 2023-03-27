library(ComplexHeatmap)
library(RColorBrewer)
data.FD = read.csv('PTRC_Master_Tumor List_20210212.csv',check.names = F)
data.FD = data.FD[data.FD$`Patient ID`!='#N/A',]

data.FZ = read.csv('PTRC_Master_Tumor List_20210212_FZ.csv',check.names = F)
data.FZ = data.FZ[!is.na(data.FZ$`Patient ID`),]

data.FP = read.csv('PTRC_Master_Tumor List_20210212_FP.csv',check.names = F)
data.FP = data.FP[!is.na(data.FP$`Patient ID`),]

sample.FD = data.frame(sample = as.character(data.FD$`Patient ID`),
                       proteome = 1,
                       phospho = 1,
                       WGS = (data.FD$`# Included in CNV Data table`>0)^2,
                       RNAseq = (data.FD$`# Included in RNAseq data table`>0)^2,
                       Tumor.Location = data.FD$`Tumor Location Group`,
                       Response = data.FD$`Tumor response`)

sample.FD$Tumor.Location = substr(sample.FD$Tumor.Location,1,2)
sample.FD$Tumor.Location[!sample.FD$Tumor.Location%in%c('OV','OM','GI')] = 'Others'

sample.FP = data.frame(sample = as.character(data.FP$`Patient ID`),
                       proteome = 1,
                       phospho = 1,
                       WGS = (data.FP$`Included in FP CNV Data table`>0)^2,
                       RNAseq = (data.FP$`Included in FP RNAseq Data table`>0)^2,
                       Tumor.Location = data.FP$`Tumor Location Group`,
                       Response = data.FP$`Tumor response`)

sample.FP$Tumor.Location = substr(sample.FP$Tumor.Location,1,2)
sample.FP$Tumor.Location[!sample.FP$Tumor.Location%in%c('OV','OM','GI')] = 'Others'

sample.FZ = data.frame(sample = as.character(data.FZ$`Patient ID`),
                       proteome = 1,
                       phospho = 1,
                       WGS = 0,
                       RNAseq = 0,
                       Tumor.Location = data.FZ$`Tumor Location Group`,
                       Response = data.FZ$`Tumor response`)

sample.FZ$Tumor.Location = substr(sample.FZ$Tumor.Location,1,2)
sample.FZ$Tumor.Location[!sample.FZ$Tumor.Location%in%c('OV','OM','GI')] = 'Others'

sample.FZ$Response = tolower(sample.FZ$Response)
sample.FD$Response = tolower(sample.FD$Response)
sample.FP$Response = tolower(sample.FP$Response)

sample.FZ = sample.FZ[order(as.numeric(as.factor(sample.FZ$Response))*1000+
                              as.numeric(as.factor(gsub('Others','A',sample.FZ$Tumor.Location)))*100+
                              sample.FZ$WGS*1+sample.FZ$RNAseq*10,decreasing = T),]

sample.FD = sample.FD[order(as.numeric(as.factor(sample.FD$Response))*1000+
                              as.numeric(as.factor(gsub('Others','A',sample.FD$Tumor.Location)))*100+
                              sample.FD$WGS*1+sample.FD$RNAseq*10,decreasing = T),]

sample.FP = sample.FP[order(as.numeric(as.factor(sample.FP$Response))*1000+
                              5-as.numeric(as.factor(gsub('Others','A',sample.FP$Tumor.Location)))*100+
                              sample.FP$WGS*1+sample.FP$RNAseq*10,decreasing = T),]

sample.all = rbind(sample.FD,sample.FP,sample.FZ)

sample.all$phospho = sample.all$phospho*2
sample.all$WGS = sample.all$WGS*3
sample.all$RNAseq = sample.all$RNAseq*4

color_sample = c('white',brewer.pal(8,'BrBG')[c(1,2,7,8)],'gold2','royalblue4')
names(color_sample) = c(0:4,'sensitive','refractory')

ha.sample = HeatmapAnnotation(Response = sample.all$Response,
                              col = list(Response = color_sample),
                              annotation_legend_param = list(direction = "horizontal",nrow = 1,title_position = "leftcenter"),
                              gp = gpar(col = "grey"))


hm.sample = Heatmap(t(as.matrix(sample.all[,2:5])),col = color_sample[1:5],
                    top_annotation = ha.sample,
                    show_column_names = F,
                    column_split = rep(c("FFPE Discovery","FFPE Validation",'Frozen Validation'), 
                                       c(dim(sample.FD)[1],dim(sample.FP)[1],dim(sample.FZ)[1])),
                    cluster_rows = F,cluster_columns = F,show_heatmap_legend = F,border = TRUE,
                    column_gap = unit(5, "mm"),
                    column_title_side = 'bottom',
                    rect_gp = gpar(type = "none"),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.rect(x = x, y = y, width = width, height = height, 
                                gp = gpar(col = "grey", fill = color_sample[as.character(t(as.matrix(sample.all[,2:5]))[i,j])]))          
                    }
)


pdf('Fig1A.pdf',height = 2.2,width = 12)
draw(hm.sample,annotation_legend_side = "bottom")
dev.off()


