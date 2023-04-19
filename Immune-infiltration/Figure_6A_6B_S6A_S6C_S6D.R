############################################################
## title: "Figure 6A_6B_S6A_S6C_S6D_S6E_S6F"
## author: "Oscar Murillo"
## date: "March 22, 2023"
## panels: "6A, 6B, S6A, S6C, S6D, S6E, S6F"
############################################################
############################################################
## Stage 0: Read Stage 0 Results
############################################################
GSE154600 = read.table("input_scRNA/GSE154600_pseudo_transformed.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
GSE154600.Markers = c("EPCAM","DCN","AIF1","EGFR","KRT14","ITGA6","KRT5","TP63","KRT17","MME","KRT8","KRT18","KRT19","FOXA1","GATA3","MUC1","CD24","KIT","GABRP","FAP","COL1A1","COL3A1","COL5A1","ACTA2","TAGLN","LUM","FBLN1","COL6A3","COL1A2","COL6A1","COL6A2","PECAM1","VWF","CDH5","SELE","PTPRC","CD2","CD3D","CD3E","CD3G","CD8A","CD8B","CD79A","CD79B","BLNK","CD14","CD68","CD163","CSF1R")
GSE154600.Subtypes = gsub(pattern = "\\.P\\d+", replacement = "", colnames(GSE154600))
GSE154600.Colors = rep("white", length(GSE154600.Subtypes))
GSE154600.Colors[GSE154600.Subtypes=="ADI"] = "red4"
GSE154600.Colors[GSE154600.Subtypes=="BCL"] = "purple"
GSE154600.Colors[GSE154600.Subtypes=="END"] = "green"
GSE154600.Colors[GSE154600.Subtypes=="EPI"] = "blue"
GSE154600.Colors[GSE154600.Subtypes=="FIB"] = "red"
GSE154600.Colors[GSE154600.Subtypes=="MAC"] = "cyan4"
GSE154600.Colors[GSE154600.Subtypes=="MON"] = "cyan"
GSE154600.Colors[GSE154600.Subtypes=="NKC"] = "yellow3"
GSE154600.Colors[GSE154600.Subtypes=="TC8"] = "orange4"

chosenGenes = read.table("input/chosen_genes.txt", header = FALSE, check.names = FALSE, sep = "\t", row.names = NULL)
chosenGenes = chosenGenes$V1

############################################################
## Read Ovarian Cancer Datasets
############################################################
FD.protein = read.table("input_protein/FD_GLBLprot_MI_FDbridge_Abund_20201002_Imput_pre_Ave.tsv", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = 1)
FD.protein = 2^FD.protein
FZ.protein = read.table("input_protein/FZ_GLBLprot_MI_FZbridge_Abund_20201002_Imput_AlignedToFD_20210823.tsv", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = 1)
FZ.protein = 2^FZ.protein

############################################################
## Read in Metadata of Datasets
############################################################
FD.meta = read.table("input_metadata/FD_GLBLprot_metadata.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)
FZ.meta = read.table("input_metadata/FZ_GLBLprot_metadata.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)

## Define FD Metadata
FD.meta.Type = FD.meta$Type
FD.meta.Response = FD.meta$Response
FD.meta.Location = FD.meta$Location
FD.meta.Neo = FD.meta$Neo
FD.meta.Source = FD.meta$Source
FD.meta.Cluster = FD.meta$Cluster

## Define FZ Metadata
FZ.meta.Type = FZ.meta$Type
FZ.meta.Response = FZ.meta$Response
FZ.meta.Location = FZ.meta$Location
FZ.meta.Cluster = FZ.meta$Cluster

############################################################
## Transform FD and FZ Datasets
############################################################
## Transform FD
input.FD.protein = FD.protein
input.FD.protein.logistic = t(apply(t(input.FD.protein), 2, logistic.a100))
input.FD.protein.logistic.nonZero =  input.FD.protein.logistic[rowSums(input.FD.protein.logistic[,]) > 0, ]
input.FD.protein.logistic.nonZero.max = (1/max(input.FD.protein.logistic.nonZero)) * input.FD.protein.logistic.nonZero

FD.protein.trans = input.FD.protein.logistic.nonZero.max
FD.protein.trans.DF = as.data.frame(FD.protein.trans)
FD.protein.trans.DF$average = rowMeans(FD.protein.trans.DF)
FD.Expressed = FD.protein.trans.DF[,"average", drop = FALSE]
FD.Expressed = FD.Expressed[FD.Expressed$average >= 0.01,,drop = FALSE]
FD.Expressed = rownames(FD.Expressed)
FD.protein.trans = FD.protein.trans[FD.Expressed,]

## Transform FZ
input.FZ.protein = FZ.protein
input.FZ.protein.logistic = t(apply(t(input.FZ.protein), 2, logistic.a100))
input.FZ.protein.logistic.nonZero =  input.FZ.protein.logistic[rowSums(input.FZ.protein.logistic[,]) > 0, ]
input.FZ.protein.logistic.nonZero.max = (1/max(input.FZ.protein.logistic.nonZero)) * input.FZ.protein.logistic.nonZero

FZ.protein.trans = input.FZ.protein.logistic.nonZero.max
FZ.protein.trans.DF = as.data.frame(FZ.protein.trans)
FZ.protein.trans.DF$average = rowMeans(FZ.protein.trans.DF)
FZ.Expressed = FZ.protein.trans.DF[,"average", drop = FALSE]
FZ.Expressed = FZ.Expressed[FZ.Expressed$average >= 0.01,,drop = FALSE]
FZ.Expressed = rownames(FZ.Expressed)
FZ.protein.trans = FZ.protein.trans[FZ.Expressed,]

## Intersect cell type specific genes with protein coverage
chosenGenes.FD = intersect(chosenGenes, rownames(FD.protein.trans))
chosenGenes.FZ = intersect(chosenGenes, rownames(FZ.protein.trans))

GSE154600.chosenGenes.FD = GSE154600[chosenGenes.FD,]
GSE154600.chosenGenes.FZ = GSE154600[chosenGenes.FZ,]

png("Stage_0/GSE154600_chosenGenes_FD.png", 1000, 1000)
heatmap.2(as.matrix(GSE154600.chosenGenes.FD),
          trace = "none", col = my_palette_1, breaks = col_breaks_1, key = FALSE, labRow = FALSE,
          ColSideColors = GSE154600.Colors, labCol = colnames(GSE154600.chosenGenes.FD), margins = c(10,10))
dev.off()

png("Stage_0/GSE154600_chosenGenes_FZ.png", 1000, 1000)
heatmap.2(as.matrix(GSE154600.chosenGenes.FZ),
          trace = "none", col = my_palette_1, breaks = col_breaks_1, key = FALSE, labRow = FALSE,
          ColSideColors = GSE154600.Colors, labCol = colnames(GSE154600.chosenGenes.FZ), margins = c(10,10))
dev.off()

############################################################
## Stage 1: Stability
############################################################
k.cells = c(3:10)
Stability.FD = estimate_stability(FD.protein.trans, chosenGenes.FD, k.cells, subset_prop = 0.8, num_subsets = 3, reps_per_subset = 2, max_its = 1000, rss_diff_stop = 1e-08)
Stability.FZ = estimate_stability(FZ.protein.trans, chosenGenes.FZ, k.cells, subset_prop = 0.8, num_subsets = 3, reps_per_subset = 2, max_its = 1000, rss_diff_stop = 1e-08)

############################################################
## Stage 1: FFPE Discovery (FD)
############################################################
FD.cellnum = Stability.FD$most_stable_num_ct

XDec.FD = run_edec_stage_1(FD.protein.trans, chosenGenes.FD, FD.cellnum, max_its = 2000, rss_diff_stop = 1e-10)
XDec.FD.Pro = XDec.FD$methylation
XDec.FD.Pro.chosenGenes = XDec.FD.Pro[rownames(XDec.FD.Pro) %in% chosenGenes.FD,]
Prop.FD = as.data.frame(XDec.FD$proportions)
Prop.FD[Prop.FD < 0] = 0

XDec.FD.names = NULL
for(i in 1:FD.cellnum){XDec.FD.names = c(XDec.FD.names, paste("FD_Profile_",i,sep = ""))}
colnames(XDec.FD.Pro.chosenGenes) = XDec.FD.names
colnames(Prop.FD) = XDec.FD.names
GSE154600.chosenGenes.FD = GSE154600[rownames(XDec.FD.Pro.chosenGenes),]
rownames(GSE154600.chosenGenes.FD) == rownames(XDec.FD.Pro.chosenGenes)

Corr.FD = cor(GSE154600.chosenGenes.FD, XDec.FD.Pro.chosenGenes, use = "complete", method = "spearman")
png("Stage_1/FD_Corr.png",1000,1000)
heatmap.2(as.matrix(Corr.FD),
          trace = "none", key = FALSE, col = my_palette_1, breaks = col_breaks_1,
          RowSideColors = GSE154600.Colors, labCol = colnames(Corr.FD), cexCol = 2, margins = c(10,10))
dev.off()

## Prepare Boxplots for FD
Prop.FD.Box = as.data.frame(Prop.FD)
colnames(Prop.FD.Box) = rownames(Corr.FD)[apply(Corr.FD, 2, which.max)]
colnames(Prop.FD.Box) = c("Epithelial","Stroma.2","Stroma.1","Immune")
Prop.FD.Box$All_Stroma = Prop.FD.Box$Stroma.1 + Prop.FD.Box$Stroma.2
Prop.FD.Box$Sample = rownames(Prop.FD.Box)
Prop.FD.Box$Type = FD.meta.Type
Prop.FD.Box$Response = FD.meta.Response
Prop.FD.Box$Source = FD.meta.Source
Prop.FD.Box$Neo = FD.meta.Neo
Prop.FD.Box$Cluster = FD.meta.Cluster
Prop.FD.Box$Location = FD.meta.Location
Prop.FD.Box$Location = gsub("Mix","Other",Prop.FD.Box$Location)
Prop.FD.Box$Location = gsub("GIOM","Other",Prop.FD.Box$Location)
Prop.FD.Box$Location = gsub("UT","Other",Prop.FD.Box$Location)
Prop.FD.Box$Location = gsub("PT","Other",Prop.FD.Box$Location)
Prop.FD.Box$Location = gsub("LN","Other",Prop.FD.Box$Location)
Prop.FD.Box$Location = gsub("GI","Other",Prop.FD.Box$Location)

Prop.FD.Box.df <- Prop.FD.Box %>% dplyr::select(Sample,Epithelial,Immune,Stroma.1,Stroma.2, everything())

############################################################
## Proportions used in Figure 6A
############################################################
write.table(Prop.FD.Box, file = "Stage_1/FD_Proportions_20210928.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Figure S6C
############################################################
Prop.FD.Box.df.Figure = Prop.FD.Box.df
Prop.FD.Box.df.Figure$All_Stroma = NULL
Prop.FD.Box.df.Figure.melt = melt(Prop.FD.Box.df.Figure)
Prop.FD.Box.df.Figure.melt$variable = gsub("Stroma.1","Stroma 1 Factor",Prop.FD.Box.df.Figure.melt$variable)
Prop.FD.Box.df.Figure.melt$variable = gsub("Stroma.2","Stroma 2 Factor",Prop.FD.Box.df.Figure.melt$variable)
Prop.FD.Box.df.Figure.melt$variable = gsub("Immune","Immune Factor",Prop.FD.Box.df.Figure.melt$variable)
Prop.FD.Box.df.Figure.melt$variable = gsub("Epithelial","Epithelial Factor",Prop.FD.Box.df.Figure.melt$variable)
Prop.FD.Box.df.Figure.melt$Cluster = gsub("Cluster_","",Prop.FD.Box.df.Figure.melt$Cluster)
Prop.FD.Box.df.Figure.melt$Cluster = factor(Prop.FD.Box.df.Figure.melt$Cluster, levels = c(1,2,3,4,5))

ggboxplot(Prop.FD.Box.df.Figure.melt, x = "Cluster", y = "value", fill = "Response", ylab = "proportions", xlab = "subtype",
          ylim = c(0,1), palette = "jco", facet.by = "variable", short.panel.labs = TRUE) +
  facet_wrap(~variable, ncol = 4) +
  theme(text = element_text(size = 16)) + 
  stat_compare_means(aes(group = Response), method = "t.test", label = "p.format", label.y = c(0.9), hide.ns = TRUE, size = 6,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))

############################################################
## Figure S6A
############################################################
## Compare to H&E Estimates
FD.HandE.Pathologist = read.table("input_HandE/Pathologist.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)

FD.HandE.Pathologist.Tumor.Ave = FD.HandE.Pathologist %>% 
  group_by(Sample) %>%
  summarise(Epithelium_Cells = list(mean(Epithelium_Cells))) %>% unnest(cols = Epithelium_Cells)
FD.HandE.Pathologist.Stroma.Ave = FD.HandE.Pathologist %>% 
  group_by(Sample) %>%
  summarise(Stroma_Cells = list(mean(Stroma_Cells))) %>% unnest(cols = Stroma_Cells)
FD.HandE.Pathologist.Lymph.Ave = FD.HandE.Pathologist %>% 
  group_by(Sample) %>%
  summarise(inflammatory_Cells = list(mean(inflammatory_Cells))) %>% unnest(cols = inflammatory_Cells)
FD.HandE.Pathologist.Adipose.Ave = FD.HandE.Pathologist %>% 
  group_by(Sample) %>%
  summarise(Adipose_Cells = list(mean(Adipose_Cells))) %>% unnest(cols = Adipose_Cells)

FD.Ave.Pathologist = merge(FD.HandE.Pathologist.Tumor.Ave, FD.HandE.Pathologist.Stroma.Ave, by = "Sample")
FD.Ave.Pathologist = merge(FD.Ave.Pathologist, FD.HandE.Pathologist.Lymph.Ave, by = "Sample")
FD.Ave.Pathologist = merge(FD.Ave.Pathologist, FD.HandE.Pathologist.Adipose.Ave, by = "Sample")
rownames(FD.Ave.Pathologist) = FD.Ave.All$Sample

Prop.FD.Box.HandE.Pathologist = merge(Prop.FD.Box, FD.Ave.Pathologist, by = "Sample")

FD.Pathologist.Epi = ggplot(Prop.FD.Box.HandE.Pathologist, aes(x = Epithelial, y = Epithelium_Cells)) +
  geom_point(size = 3) + xlab("Epithelial %") + ylab("Pathologist epithelium cells %") +
  geom_smooth(method = lm, fullrange = FALSE) + stat_cor(label.y = 1.25, size = 6) +
  theme_bw() + theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

FD.Pathologist.Stroma = ggplot(Prop.FD.Box.HandE.Pathologist, aes(x = All_Stroma, y = Stroma_Cells+Adipose_Cells)) +
  geom_point(size = 3) + xlab("Stroma %") + ylab("Pathologist stromal cells %") +
  geom_smooth(method = lm, fullrange = FALSE) + stat_cor(size = 6) + 
  theme_bw() + theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

FD.Pathologist.Immune = ggplot(Prop.FD.Box.HandE.Pathologist, aes(x = Immune, y = inflammatory_Cells)) +
  geom_point(size = 3) + xlab("Immune %") + ylab("Pathologist inflammatory cells %") +
  geom_smooth(method = lm, fullrange = FALSE) + stat_cor(size = 6) +
  theme_bw() + theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggarrange(FD.Pathologist.Epi, FD.Pathologist.Stroma, FD.Pathologist.Immune, ncol = 3, nrow = 1)

############################################################
## Stage 1: Frozen Discovery (FZ)
############################################################
FZ.cellnum = Stability.FZ$most_stable_num_ct

XDec.FZ = run_edec_stage_1(FZ.protein.trans, chosenGenes.FZ, FZ.cellnum, max_its = 2000, rss_diff_stop = 1e-10)
XDec.FZ.Pro = XDec.FZ$methylation
XDec.FZ.Pro.chosenGenes = XDec.FZ.Pro[rownames(XDec.FZ.Pro) %in% chosenGenes.FZ,]
Prop.FZ = as.data.frame(XDec.FZ$proportions)
Prop.FZ[Prop.FZ < 0] = 0

XDec.FZ.names = NULL
for(i in 1:FZ.cellnum){XDec.FZ.names = c(XDec.FZ.names, paste("FZ_Profile_",i,sep = ""))}
colnames(XDec.FZ.Pro.chosenGenes) = XDec.FZ.names
colnames(Prop.FZ) = XDec.FZ.names
GSE154600.chosenGenes.FZ = GSE154600[rownames(XDec.FZ.Pro.chosenGenes),]
rownames(GSE154600.chosenGenes.FZ) == rownames(XDec.FZ.Pro.chosenGenes)

Corr.FZ = cor(GSE154600.chosenGenes.FZ, XDec.FZ.Pro.chosenGenes, use = "complete", method = "spearman")
heatmap.2(as.matrix(Corr.FZ),
          trace = "none", key = FALSE, col = my_palette_1, breaks = col_breaks_1,
          RowSideColors = GSE154600.Colors, labCol = colnames(Corr.FZ), cexCol = 2, margins = c(10,10))

## Prepare Boxplots for FZ
Prop.FZ.Box = as.data.frame(Prop.FZ)
colnames(Prop.FZ.Box) = rownames(Corr.FZ)[apply(Corr.FZ, 2, which.max)]
colnames(Prop.FZ.Box) = c("Stroma.2","Stroma.1","Epithelial","Immune")
Prop.FZ.Box$All_Stroma = Prop.FZ.Box$Stroma.1 + Prop.FZ.Box$Stroma.2
Prop.FZ.Box$Sample = rownames(Prop.FZ.Box)
Prop.FZ.Box$Type = FZ.meta.Type
Prop.FZ.Box$Response = FZ.meta.Response
Prop.FZ.Box$Cluster = FZ.meta.Cluster
Prop.FZ.Box$Location = FZ.meta.Location
Prop.FZ.Box$Location = gsub("PT","Other",Prop.FZ.Box$Location)
Prop.FZ.Box$Location = gsub("LN","Other",Prop.FZ.Box$Location)
Prop.FZ.Box$Location = gsub("GI","Other",Prop.FZ.Box$Location)

Prop.FZ.Box.df <- Prop.FZ.Box %>% dplyr::select(Sample,Epithelial,Immune,Stroma.1,Stroma.2, everything())
write.table(Prop.FZ.Box.df, file = "Stage_1/FZ_Proportions_20210928.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## FZ.pathway.predict.order is from previous analysis estimating PAM clusters for FZ
## Plot Proportions across Clusters
FZ.meta.CombinePathway = merge(Prop.FZ.Box.df, FZ.pathway.predict.order[,c(1,7)], by = "Sample")
colnames(FZ.meta.CombinePathway)[11] = c("PAM_Pathway")
Prop.FZ.Box.df.Figure = FZ.meta.CombinePathway
Prop.FZ.Box.df.Figure$All_Stroma = NULL
Prop.FZ.Box.df.Figure.melt = melt(Prop.FZ.Box.df.Figure)
Prop.FZ.Box.df.Figure.melt$variable = gsub("Stroma.1","Stroma 1 Factor",Prop.FZ.Box.df.Figure.melt$variable)
Prop.FZ.Box.df.Figure.melt$variable = gsub("Stroma.2","Stroma 2 Factor",Prop.FZ.Box.df.Figure.melt$variable)
Prop.FZ.Box.df.Figure.melt$variable = gsub("Immune","Immune Factor",Prop.FZ.Box.df.Figure.melt$variable)
Prop.FZ.Box.df.Figure.melt$variable = gsub("Epithelial","Epithelial Factor",Prop.FZ.Box.df.Figure.melt$variable)
Prop.FZ.Box.df.Figure.melt$PAM_Pathway = gsub("Cluster_","",Prop.FZ.Box.df.Figure.melt$PAM_Pathway)
Prop.FZ.Box.df.Figure.melt$PAM_Pathway = factor(Prop.FZ.Box.df.Figure.melt$PAM_Pathway, levels = c(1,2,3,4,5))
Prop.FZ.Box.df.Figure.melt = Prop.FZ.Box.df.Figure.melt[Prop.FZ.Box.df.Figure.melt$PAM_Pathway %in% c(1,3,5),]

############################################################
## Figure S6D
############################################################
ggplot(Prop.FZ.Box.df.Figure.melt, aes(x = PAM_Pathway, y = value, color = Response)) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge()) +
  facet_wrap(~variable, ncol = 4) + theme_classic() +
  labs(x = "subtype", y = "proportion") + scale_color_manual(values = c("#0073C2","#EFC000")) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 0.5, hide.ns = TRUE) + scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

############################################################
## Figure 6B
############################################################
Prop.FZ.Box.df.Figure.melt.All = Prop.FZ.Box.df.Figure.melt
Prop.FZ.Box.df.Figure.melt.All = Prop.FZ.Box.df.Figure.melt.All[Prop.FZ.Box.df.Figure.melt.All$variable == "Immune Factor",]
Prop.FZ.Box.df.Figure.melt.All = Prop.FZ.Box.df.Figure.melt.All[Prop.FZ.Box.df.Figure.melt.All$PAM_Pathway == 5,]
Prop.FZ.Box.df.Figure.melt.All = Prop.FZ.Box.df.Figure.melt.All[,c(1,3,6,8)]
colnames(Prop.FZ.Box.df.Figure.melt.All) = c("Sample","Response","subtype","immune")
Prop.FZ.Box.df.Figure.melt.All$cohort = c("FZ All")

Prop.FZ.Box.df.Figure.melt.NonOverlap = Prop.FZ.Box.df.Figure.melt[!Prop.FZ.Box.df.Figure.melt$Sample %in% FZ.sampleOverlap,]
Prop.FZ.Box.df.Figure.melt.NonOverlap = Prop.FZ.Box.df.Figure.melt.NonOverlap[Prop.FZ.Box.df.Figure.melt.NonOverlap$variable == "Immune Factor",]
Prop.FZ.Box.df.Figure.melt.NonOverlap = Prop.FZ.Box.df.Figure.melt.NonOverlap[Prop.FZ.Box.df.Figure.melt.NonOverlap$PAM_Pathway == 5,]
Prop.FZ.Box.df.Figure.melt.NonOverlap = Prop.FZ.Box.df.Figure.melt.NonOverlap[,c(1,3,6,8)]
colnames(Prop.FZ.Box.df.Figure.melt.NonOverlap) = c("Sample","Response","subtype","immune")
Prop.FZ.Box.df.Figure.melt.NonOverlap$cohort = c("FZ NonOverlap")

Prop.FZ.Box.df.Figure.melt.Overlap = Prop.FZ.Box.df.Figure.melt[Prop.FZ.Box.df.Figure.melt$Sample %in% FZ.sampleOverlap,]
Prop.FZ.Box.df.Figure.melt.Overlap = Prop.FZ.Box.df.Figure.melt.Overlap[Prop.FZ.Box.df.Figure.melt.Overlap$variable == "Immune Factor",]
Prop.FZ.Box.df.Figure.melt.Overlap = Prop.FZ.Box.df.Figure.melt.Overlap[Prop.FZ.Box.df.Figure.melt.Overlap$PAM_Pathway == 5,]
Prop.FZ.Box.df.Figure.melt.Overlap = Prop.FZ.Box.df.Figure.melt.Overlap[,c(1,3,6,8)]
colnames(Prop.FZ.Box.df.Figure.melt.Overlap) = c("Sample","Response","subtype","immune")
Prop.FZ.Box.df.Figure.melt.Overlap$cohort = c("FZ Overlap")

Prop.FD.Box.Immune = Prop.FD.Box[,c("Immune","Sample","Response","Cluster")]
Prop.FD.Box.Immune = Prop.FD.Box.Immune[Prop.FD.Box.Immune$Cluster == "Cluster_5",]
Prop.FD.Box.Immune.Figure = Prop.FD.Box.Immune[,c(2,3,4,1)]
colnames(Prop.FD.Box.Immune.Figure) = c("Sample","Response","subtype","immune")
Prop.FD.Box.Immune.Figure$subtype = gsub("Cluster_","",Prop.FD.Box.Immune.Figure$subtype)
rownames(Prop.FD.Box.Immune.Figure) = NULL
Prop.FD.Box.Immune.Figure$cohort = c("FFPE")

Prop.Immune.Figure = rbind(Prop.FD.Box.Immune.Figure, Prop.FZ.Box.df.Figure.melt.Overlap, Prop.FZ.Box.df.Figure.melt.NonOverlap, Prop.FZ.Box.df.Figure.melt.All)
Prop.Immune.Figure.melt = melt(Prop.Immune.Figure)
Prop.Immune.Figure.melt$cohort = factor(Prop.Immune.Figure.melt$cohort, levels = c("FFPE","FZ Overlap","FZ NonOverlap","FZ All"))
Prop.Immune.Figure.melt = Prop.Immune.Figure.melt[Prop.Immune.Figure.melt$cohort %in% c("FFPE","FZ NonOverlap"),]

ggplot(Prop.Immune.Figure.melt, aes(x = variable, y = value, color = Response)) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge()) +
  facet_wrap(~cohort, scales = "free_y", ncol = 4) + scale_color_manual(values = c("#0073C2","#E7B800")) +
  labs(title = "subtype 5", x = "", y = "proportion (%)") + theme_classic() +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 0.5, vjust = -1) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))
