############################################################
## title: "Figure 7C_7D_S7D_S7G"
## author: "Oscar Murillo"
## date: "March 22, 2023"
## panels: "7C, 7D, S7D, S7G"
############################################################
############################################################
## FD - Peptide Analysis
############################################################
## Read in All Peptide Search 140067
FD.AllPeptide.Test = read.table("peptide_analysis/PTRC_FD_noKR_peptide_ratio.txt", header = TRUE, quote = "", check.names = FALSE, sep = "\t", row.names = NULL)
FD.AllPeptide.Test = FD.AllPeptide.Test[rowSums(is.na(FD.AllPeptide.Test)[,14:244]) != 231,]
FD.AllPeptide.Test = FD.AllPeptide.Test[,c(14:244)]
FD.AllPeptide.Test.med = unlist(lapply(FD.AllPeptide.Test, median, na.rm = TRUE))
names(FD.AllPeptide.Test.med) = FD.PlexSamples
FD.AllPeptide.Test.med = FD.AllPeptide.Test.med[-grep("Cell", names(FD.AllPeptide.Test.med))]
FD.AllPeptide.Test.med = FD.AllPeptide.Test.med[-grep("Frozen", names(FD.AllPeptide.Test.med))]
FD.AllPeptide.Test.med = FD.AllPeptide.Test.med[-grep("FFPE", names(FD.AllPeptide.Test.med))]

FD.AllPeptide.Test.med["FD_p01_UA_CU-23_S-15-0033313_B1"] = (FD.AllPeptide.Test.med["FD_p01_UA_CU-23_S-15-0033313_B1"] + FD.AllPeptide.Test.med["FD_p09_UA_CU-23_S-15-0033313_B1"] + FD.AllPeptide.Test.med["FD_p17_UA_CU-23_S-15-0033313_B1"])/3
FD.AllPeptide.Test.med["FD_p02_MC_7598_40885_B1"] = (FD.AllPeptide.Test.med["FD_p02_MC_7598_40885_B1"] + FD.AllPeptide.Test.med["FD_p10_MC_7598_40885_B1"] + FD.AllPeptide.Test.med["FD_p18_MC_7598_40885_B1"])/3
FD.AllPeptide.Test.med["FD_p04_UA_CU-15_S-13-0034287_A2"] = (FD.AllPeptide.Test.med["FD_p04_UA_CU-15_S-13-0034287_A2"] + FD.AllPeptide.Test.med["FD_p12_UA_CU-15_S-13-0034287_A2"] + FD.AllPeptide.Test.med["FD_p19_UA_CU-15_S-13-0034287_A2"])/3
FD.AllPeptide.Test.med["FD_p06_FH_94_654490"] = (FD.AllPeptide.Test.med["FD_p06_FH_94_654490"] + FD.AllPeptide.Test.med["FD_p14_FH_94_654490"] + FD.AllPeptide.Test.med["FD_p20_FH_94_654490"])/3
FD.AllPeptide.Test.med["FD_p08_UA_CU-22_S-15-0030998_A1"] = (FD.AllPeptide.Test.med["FD_p08_UA_CU-22_S-15-0030998_A1"] + FD.AllPeptide.Test.med["FD_p16_UA_CU-22_S-15-0030998_A1"] + FD.AllPeptide.Test.med["FD_p21_UA_CU-22_S-15-0030998_A1"])/3
FD.AllPeptide.Test.med.Test = as.data.frame(FD.AllPeptide.Test.med)
FD.AllPeptide.Test.med.Test$Sample = rownames(FD.AllPeptide.Test.med.Test)
FD.AllPeptide.Test.med.Test = FD.AllPeptide.Test.med.Test[!FD.AllPeptide.Test.med.Test$Sample %in% FD.NonKR.Test.Remove,]
FD.AllPeptide.Test.med.Median = FD.AllPeptide.Test.med.Test$FD.AllPeptide.Test.med

FD.NonKR.Test = FD.AllPeptide.Search.Non
rownames(FD.NonKR.Test) = FD.NonKR.Test$Sequence
FD.NonKR.Test = FD.NonKR.Test[,c(14:244)]
colnames(FD.NonKR.Test) = FD.PlexSamples
FD.NonKR.Test = FD.NonKR.Test[-grep("Cell", names(FD.NonKR.Test))]
FD.NonKR.Test = FD.NonKR.Test[-grep("Frozen", names(FD.NonKR.Test))]
FD.NonKR.Test = FD.NonKR.Test[-grep("FFPE", names(FD.NonKR.Test))]

FD.NonKR.Test["FD_p01_UA_CU-23_S-15-0033313_B1"] = (FD.NonKR.Test["FD_p01_UA_CU-23_S-15-0033313_B1"] + FD.NonKR.Test["FD_p09_UA_CU-23_S-15-0033313_B1"] + FD.NonKR.Test["FD_p17_UA_CU-23_S-15-0033313_B1"])/3
FD.NonKR.Test["FD_p02_MC_7598_40885_B1"] = (FD.NonKR.Test["FD_p02_MC_7598_40885_B1"] + FD.NonKR.Test["FD_p10_MC_7598_40885_B1"] + FD.NonKR.Test["FD_p18_MC_7598_40885_B1"])/3
FD.NonKR.Test["FD_p04_UA_CU-15_S-13-0034287_A2"] = (FD.NonKR.Test["FD_p04_UA_CU-15_S-13-0034287_A2"] + FD.NonKR.Test["FD_p12_UA_CU-15_S-13-0034287_A2"] + FD.NonKR.Test["FD_p19_UA_CU-15_S-13-0034287_A2"])/3
FD.NonKR.Test["FD_p06_FH_94_654490"] = (FD.NonKR.Test["FD_p06_FH_94_654490"] + FD.NonKR.Test["FD_p14_FH_94_654490"] + FD.NonKR.Test["FD_p20_FH_94_654490"])/3
FD.NonKR.Test["FD_p08_UA_CU-22_S-15-0030998_A1"] = (FD.NonKR.Test["FD_p08_UA_CU-22_S-15-0030998_A1"] + FD.NonKR.Test["FD_p16_UA_CU-22_S-15-0030998_A1"] + FD.NonKR.Test["FD_p21_UA_CU-22_S-15-0030998_A1"])/3
FD.NonKR.Test = as.data.frame(FD.NonKR.Test)

FD.NonKR.Test.Remove = c("FD_p09_UA_CU-23_S-15-0033313_B1","FD_p17_UA_CU-23_S-15-0033313_B1","FD_p10_MC_7598_40885_B1","FD_p18_MC_7598_40885_B1","FD_p12_UA_CU-15_S-13-0034287_A2","FD_p19_UA_CU-15_S-13-0034287_A2","FD_p14_FH_94_654490","FD_p20_FH_94_654490","FD_p16_UA_CU-22_S-15-0030998_A1","FD_p21_UA_CU-22_S-15-0030998_A1")
FD.NonKR.Test = FD.NonKR.Test[,!colnames(FD.NonKR.Test) %in% FD.NonKR.Test.Remove]
FD.NonKR.Test.Norm = sweep(FD.NonKR.Test, 2, FD.AllPeptide.Test.med.Median, `/`)

FD.NonKR.Test.NA = FD.NonKR.Test.Norm[rowSums(is.na(FD.NonKR.Test.Norm)) < 118,]
FD.NonKR.Test.NA[is.na(FD.NonKR.Test.NA)] <- 1

FD.NonKR.Meta = FD.Figure3.meta.order
FD.NonKR.Meta = FD.NonKR.Meta[colnames(FD.NonKR.Test),]
FD.NonKR.Meta$TMT = substr(FD.NonKR.Meta$Sample, 1, 6)

FD.NonKR.Meta.Ann = HeatmapAnnotation(site = FD.NonKR.Meta$site, sample.age = FD.NonKR.Meta$age, location = FD.NonKR.Meta$location,
                                      neo.adjuvant = FD.NonKR.Meta$neo.adjuvant, TCGA.subtype = FD.NonKR.Meta$TCGA.subtype, tumor.purity = FD.NonKR.Meta$tumor.purity,
                                      genome.instability = FD.NonKR.Meta$genome.instability, prediction.score = FD.NonKR.Meta$prediction.score,
                                      response = FD.NonKR.Meta$response, subtype = FD.NonKR.Meta$subtype, TMT = FD.NonKR.Meta$TMT)

col_fun = colorRamp2(c(0, 1, 2), c("blue", "white", "red"))

FD.NonKR.Test.NA = as.matrix(FD.NonKR.Test.NA)
theme_set(theme_grey(base_size = 16))
pdf("FD_NonKR_Test.pdf", height = 12, width = 16, useDingbats = FALSE)
ht_list = Heatmap(FD.NonKR.Test.NA, top_annotation = FD.NonKR.Meta.Ann,
                  cluster_rows = TRUE, cluster_columns = TRUE, column_split = factor(FD.NonKR.Meta$subtype, levels = c("1","2","3","4","5")),
                  border = TRUE, na_col = "white", col = col_fun,
                  show_column_names = FALSE, show_row_names = FALSE)
draw(ht_list, heatmap_legend_side = c("right"), annotation_legend_side = "right",  padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

ht_list <- draw(ht_list)
FD.row.order = row_order(ht_list)
FD.col.order = as.vector(unlist(column_order(ht_list)))

FD.NonKR.Test.NA.Row = as.data.frame(FD.NonKR.Test.NA[FD.row.order, FD.col.order])
FD.NonKR.Test.NA.Row.Bottom = FD.NonKR.Test.NA.Row[c(460:538),]
FD.NonKR.Test.NA.Row.Bottom = as.matrix(FD.NonKR.Test.NA.Row.Bottom)

FD.AllPeptide.Search.Non.Bottom = as.data.frame(FD.NonKR.Test.Norm[rownames(FD.NonKR.Test.NA.Row.Bottom),])
FD.AllPeptide.Search.Non.Bottom.med = unlist(lapply(FD.AllPeptide.Search.Non.Bottom, median, na.rm = TRUE))

Bottom.Ratio = FD.AllPeptide.Search.Non.Bottom.med
names(Bottom.Ratio) = FD.PlexSamples
Bottom.Ratio = Bottom.Ratio[-grep("Cell", names(Bottom.Ratio))]
Bottom.Ratio = Bottom.Ratio[-grep("Frozen", names(Bottom.Ratio))]
Bottom.Ratio = Bottom.Ratio[-grep("FFPE", names(Bottom.Ratio))]

Bottom.Ratio["FD_p01_UA_CU-23_S-15-0033313_B1"] = (Bottom.Ratio["FD_p01_UA_CU-23_S-15-0033313_B1"] + Bottom.Ratio["FD_p09_UA_CU-23_S-15-0033313_B1"] + Bottom.Ratio["FD_p17_UA_CU-23_S-15-0033313_B1"])/3
Bottom.Ratio["FD_p02_MC_7598_40885_B1"] = (Bottom.Ratio["FD_p02_MC_7598_40885_B1"] + Bottom.Ratio["FD_p10_MC_7598_40885_B1"] + Bottom.Ratio["FD_p18_MC_7598_40885_B1"])/3
Bottom.Ratio["FD_p04_UA_CU-15_S-13-0034287_A2"] = (Bottom.Ratio["FD_p04_UA_CU-15_S-13-0034287_A2"] + Bottom.Ratio["FD_p12_UA_CU-15_S-13-0034287_A2"] + Bottom.Ratio["FD_p19_UA_CU-15_S-13-0034287_A2"])/3
Bottom.Ratio["FD_p06_FH_94_654490"] = (Bottom.Ratio["FD_p06_FH_94_654490"] + Bottom.Ratio["FD_p14_FH_94_654490"] + Bottom.Ratio["FD_p20_FH_94_654490"])/3
Bottom.Ratio["FD_p08_UA_CU-22_S-15-0030998_A1"] = (Bottom.Ratio["FD_p08_UA_CU-22_S-15-0030998_A1"] + Bottom.Ratio["FD_p16_UA_CU-22_S-15-0030998_A1"] + Bottom.Ratio["FD_p21_UA_CU-22_S-15-0030998_A1"])/3
Bottom.Ratio = as.data.frame(FD.AllPeptide.Search.Non.Bottom.med)
Bottom.Ratio$Sample = rownames(Bottom.Ratio)  
colnames(Bottom.Ratio) = c("Bottom.Ratio","Sample")

FD.QC.meta.merge.Bottom = merge(FD.QC.meta.merge, Bottom.Ratio, by = "Sample")
FD.QC.meta.merge.Bottom$Stroma.1 = NULL
FD.QC.meta.merge.Bottom$Stroma.2 = NULL
FD.QC.meta.merge.Bottom = FD.QC.meta.merge.Bottom[,c("Sample","Response","Type","Location","Source","Cluster","Epithelial","Immune","All_Stroma","Mox","MC","Cterm","NonOverFully.Ratio","SemiOverFully.Ratio","Bottom.Ratio")]
FD.QC.meta.merge.Bottom.melt = melt(FD.QC.meta.merge.Bottom)
FD.QC.meta.merge.Bottom.melt$Cluster = gsub("Cluster_", "", FD.QC.meta.merge.Bottom.melt$Cluster)
FD.QC.meta.merge.Bottom.melt$Cluster = as.factor(FD.QC.meta.merge.Bottom.melt$Cluster)
FD.QC.meta.merge.Bottom.melt.Figure = FD.QC.meta.merge.Bottom.melt
FD.QC.meta.merge.Bottom.melt.Figure = FD.QC.meta.merge.Bottom.melt.Figure[FD.QC.meta.merge.Bottom.melt.Figure$variable %in% c("NonOverFully.Ratio","SemiOverFully.Ratio","Bottom.Ratio"),]
FD.QC.meta.merge.Bottom.melt.Figure$variable = gsub("NonOverFully.Ratio","NonKR peptides",FD.QC.meta.merge.Bottom.melt.Figure$variable)
FD.QC.meta.merge.Bottom.melt.Figure$variable = gsub("SemiOverFully.Ratio","SemiKR peptides",FD.QC.meta.merge.Bottom.melt.Figure$variable)
FD.QC.meta.merge.Bottom.melt.Figure$variable = gsub("Bottom.Ratio","Elevated NonKR peptides",FD.QC.meta.merge.Bottom.melt.Figure$variable)
FD.QC.meta.merge.Bottom.melt.Figure$variable = factor(FD.QC.meta.merge.Bottom.melt.Figure$variable, levels = c("SemiKR peptides","NonKR peptides","Elevated NonKR peptides"))

FD.QC.meta.merge.Bottom.melt.Figure.nonKR = FD.QC.meta.merge.Bottom.melt.Figure[FD.QC.meta.merge.Bottom.melt.Figure$variable == "NonKR peptides",]
FD.QC.meta.merge.Bottom.melt.Figure.SemiKR = FD.QC.meta.merge.Bottom.melt.Figure[FD.QC.meta.merge.Bottom.melt.Figure$variable == "SemiKR peptides",]
FD.QC.meta.merge.Bottom.melt.Figure.MMPnonKR = FD.QC.meta.merge.Bottom.melt.Figure[FD.QC.meta.merge.Bottom.melt.Figure$variable == "Elevated NonKR peptides",]

write.table(FD.QC.meta.merge.Bottom.melt.Figure.nonKR, "quality/FD.QC.meta.merge.Bottom.melt.Figure.nonKR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(FD.QC.meta.merge.Bottom.melt.Figure.SemiKR, "quality/FD.QC.meta.merge.Bottom.melt.Figure.SemiKR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(FD.QC.meta.merge.Bottom.melt.Figure.MMPnonKR, "quality/FD.QC.meta.merge.Bottom.melt.Figure.MMPnonKR.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Figure 7C
############################################################
theme_set(theme_grey(base_size = 20))
pdf("FD_Boxplot_Ratios_Cluster.pdf", width = 6, height = 5)
ggplot(FD.QC.meta.merge.Bottom.melt.Figure.nonKR, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Cluster), size = 2, shape = 21, position = position_jitterdodge()) +
  theme_bw() + ylab("NonKR ratio") + xlab("clusters") + stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")

ggplot(FD.QC.meta.merge.Bottom.melt.Figure.SemiKR, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Cluster), size = 2, shape = 21, position = position_jitterdodge()) +
  theme_bw() + ylab("SemiKR ratio") + xlab("clusters") + stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")

ggplot(FD.QC.meta.merge.Bottom.melt.Figure.MMPnonKR, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Cluster), size = 2, shape = 21, position = position_jitterdodge()) +
  theme_bw() + ylab("NonKR peptides enriched \nin ECM-related proteins") + xlab("clusters") + stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")
dev.off()

############################################################
## FZ - Peptide Analysis
############################################################
## Read in All Peptide Search 140067
FZ.AllPeptide.Test = read.table("peptide_analysis/peptides_ratios_fz_noKR_location.txt", header = TRUE, quote = "", check.names = FALSE, sep = "\t", row.names = NULL)
FZ.AllPeptide.Test = FZ.AllPeptide.Test[rowSums(is.na(FZ.AllPeptide.Test)[,12:99]) != 88,]
FZ.AllPeptide.Test = FZ.AllPeptide.Test[,-grep("intensity 10", colnames(FZ.AllPeptide.Test))]
FZ.AllPeptide.Test = FZ.AllPeptide.Test[,-grep("intensity 11", colnames(FZ.AllPeptide.Test))]

FZ.AllPeptide.Test = FZ.AllPeptide.Test[,c(12:83)]
FZ.AllPeptide.Test.med = unlist(lapply(FZ.AllPeptide.Test, median, na.rm = TRUE))
names(FZ.AllPeptide.Test.med) = FZ.PlexSamples

FZ.AllPeptide.Test.med["FZ_p01_FH_35_261437"] = (FZ.AllPeptide.Test.med["FZ_p01_FH_35_261437"] + FZ.AllPeptide.Test.med["FZ_p04_FH_35_261437"] + FZ.AllPeptide.Test.med["FZ_p07_FH_35_261437"])/3
FZ.AllPeptide.Test.med["FZ_p01_FH_99_439629"] = (FZ.AllPeptide.Test.med["FZ_p01_FH_99_439629"] + FZ.AllPeptide.Test.med["FZ_p03_FH_99_439629"] + FZ.AllPeptide.Test.med["FZ_p05_FH_99_439629"])/3
FZ.AllPeptide.Test.med["FZ_p02_FH_73_497549"] = (FZ.AllPeptide.Test.med["FZ_p02_FH_73_497549"] + FZ.AllPeptide.Test.med["FZ_p05_FH_73_497549"] + FZ.AllPeptide.Test.med["FZ_p08_FH_73_497549"])/3
FZ.AllPeptide.Test.med["FZ_p03_FH_98_371244"] = (FZ.AllPeptide.Test.med["FZ_p03_FH_98_371244"] + FZ.AllPeptide.Test.med["FZ_p06_FH_98_371244"] + FZ.AllPeptide.Test.med["FZ_p08_FH_98_371244"])/3
FZ.NonKR.Test.Remove = c("FZ_p04_FH_35_261437","FZ_p07_FH_35_261437","FZ_p03_FH_99_439629","FZ_p05_FH_99_439629","FZ_p05_FH_73_497549","FZ_p08_FH_73_497549","FZ_p06_FH_98_371244","FZ_p08_FH_98_371244")

FZ.AllPeptide.Test.med.Test = as.data.frame(FZ.AllPeptide.Test.med)
FZ.AllPeptide.Test.med.Test$Sample = rownames(FZ.AllPeptide.Test.med.Test)
FZ.AllPeptide.Test.med.Test = FZ.AllPeptide.Test.med.Test[!FZ.AllPeptide.Test.med.Test$Sample %in% FZ.NonKR.Test.Remove,]
FZ.AllPeptide.Test.med.Median = FZ.AllPeptide.Test.med.Test$FZ.AllPeptide.Test.med

FZ.NonKR.Test = FZ.AllPeptide.Search.Non
rownames(FZ.NonKR.Test) = FZ.NonKR.Test$Sequence
FZ.NonKR.Test = FZ.NonKR.Test[,c(12:83)]
colnames(FZ.NonKR.Test) = FZ.PlexSamples

FZ.NonKR.Test["FZ_p01_FH_35_261437"] = (FZ.NonKR.Test["FZ_p01_FH_35_261437"] + FZ.NonKR.Test["FZ_p04_FH_35_261437"] + FZ.NonKR.Test["FZ_p07_FH_35_261437"])/3
FZ.NonKR.Test["FZ_p01_FH_99_439629"] = (FZ.NonKR.Test["FZ_p01_FH_99_439629"] + FZ.NonKR.Test["FZ_p03_FH_99_439629"] + FZ.NonKR.Test["FZ_p05_FH_99_439629"])/3
FZ.NonKR.Test["FZ_p02_FH_73_497549"] = (FZ.NonKR.Test["FZ_p02_FH_73_497549"] + FZ.NonKR.Test["FZ_p05_FH_73_497549"] + FZ.NonKR.Test["FZ_p08_FH_73_497549"])/3
FZ.NonKR.Test["FZ_p03_FH_98_371244"] = (FZ.NonKR.Test["FZ_p03_FH_98_371244"] + FZ.NonKR.Test["FZ_p06_FH_98_371244"] + FZ.NonKR.Test["FZ_p08_FH_98_371244"])/3
FZ.NonKR.Test = as.data.frame(FZ.NonKR.Test)

FZ.NonKR.Test = FZ.NonKR.Test[,!colnames(FZ.NonKR.Test) %in% FZ.NonKR.Test.Remove]
FZ.NonKR.Test.Norm = sweep(FZ.NonKR.Test, 2, FZ.AllPeptide.Test.med.Median, `/`)

FZ.NonKR.Test.NA = FZ.NonKR.Test.Norm[rowSums(is.na(FZ.NonKR.Test.Norm)) < 48,]
FZ.NonKR.Test.NA[is.na(FZ.NonKR.Test.NA)] <- 1

FZ.NonKR.Meta = FZ.Figure3.meta.order
FZ.NonKR.Meta = FZ.NonKR.Meta[colnames(FZ.NonKR.Test),]
FZ.NonKR.Meta$TMT = substr(FZ.NonKR.Meta$Sample, 1, 6)

FZ.NonKR.Meta.Ann = HeatmapAnnotation(site = FZ.NonKR.Meta$Source, sample.age = FZ.NonKR.Meta$Age, location = FZ.NonKR.Meta$Location,
                                      response = FZ.NonKR.Meta$response, subtype = FZ.NonKR.Meta$subtype, TMT = FZ.NonKR.Meta$TMT)

col_fun = colorRamp2(c(0, 1, 2), c("blue", "white", "red"))

FZ.NonKR.Test.NA = as.matrix(FZ.NonKR.Test.NA)
theme_set(theme_grey(base_size = 16))
pdf("FZ_NonKR_Test.pdf", height = 12, width = 16, useDingbats = FALSE)
FZ_ht_list = Heatmap(FZ.NonKR.Test.NA, top_annotation = FZ.NonKR.Meta.Ann,
                     cluster_rows = TRUE, cluster_columns = TRUE, column_split = factor(FZ.NonKR.Meta$subtype, levels = c("1","2","3","4","5")),
                     border = TRUE, na_col = "white", col = col_fun,
                     show_column_names = FALSE, show_row_names = FALSE)
draw(FZ_ht_list, heatmap_legend_side = c("right"), annotation_legend_side = "right",  padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

FZ_ht_list <- draw(FZ_ht_list)
FZ.row.order = row_order(FZ_ht_list)
FZ.col.order = as.vector(unlist(column_order(FZ_ht_list)))

FZ.NonKR.Test.NA.Row = as.data.frame(FZ.NonKR.Test.NA[FZ.row.order, FZ.col.order])
FZ.NonKR.Test.NA.Row.Bottom = FZ.NonKR.Test.NA.Row[c(489:609,629:692),]
FZ.NonKR.Test.NA.Row.Bottom = as.matrix(FZ.NonKR.Test.NA.Row.Bottom)

FZ.NonKR.Meta.Row = as.data.frame(FZ.NonKR.Meta[FZ.col.order,])
FZ.NonKR.Meta.Row.Ann = HeatmapAnnotation(site = FZ.NonKR.Meta.Row$Source, sample.age = FZ.NonKR.Meta.Row$Age, location = FZ.NonKR.Meta.Row$Location,
                                          response = FZ.NonKR.Meta.Row$Response, subtype = FZ.NonKR.Meta.Row$subtype, TMT = FZ.NonKR.Meta.Row$TMT,
                                          col = list(site = c("FHCRC" = "turquoise1","Mayo" = "blue","UAB" = "green"),
                                                     location = c("OV" = "yellow","OM" = "red","other" = "dodgerblue","mix" = "orange"),
                                                     response = c("sensitive" = "gold","refractory" = "blue"),
                                                     subtype = c("1" = "#F8766D","2" = "#BB9D00","3" = "#00B81F","4" = "#00A5FF","5" = "#E76BF3")), na_col = "lightgray")

theme_set(theme_grey(base_size = 16))
pdf("FZ_NonKR_Test_Bottom.pdf", height = 12, width = 16, useDingbats = FALSE)
ht_list_2 = Heatmap(as.matrix(FZ.NonKR.Test.NA.Row.Bottom), top_annotation = FZ.NonKR.Meta.Row.Ann,
                    cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(FZ.NonKR.Meta.Row$subtype, levels = c("1","2","3","4","5")),
                    border = TRUE, na_col = "white", col = col_fun,
                    show_column_names = FALSE, show_row_names = TRUE)
draw(ht_list_2, heatmap_legend_side = c("right"), annotation_legend_side = "right",  padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

FZ.AllPeptide.Search.Non.FZ.Bottom = as.data.frame(FZ.NonKR.Test.Norm[rownames(FZ.NonKR.Test.NA.Row.Bottom),])
FZ.AllPeptide.Search.Non.FZ.Bottom.med = unlist(lapply(FZ.AllPeptide.Search.Non.FZ.Bottom, median, na.rm = TRUE))

FZ.Bottom.Ratio = FZ.AllPeptide.Search.Non.FZ.Bottom.med
FZ.Bottom.Ratio = as.data.frame(FZ.AllPeptide.Search.Non.FZ.Bottom.med)
FZ.Bottom.Ratio$Sample = rownames(FZ.Bottom.Ratio)  
colnames(FZ.Bottom.Ratio) = c("Bottom.Ratio","Sample")

FZ.QC.meta.merge.Bottom = merge(FZ.QC.meta.merge, FZ.Bottom.Ratio, by = "Sample")
FZ.QC.meta.merge.Bottom.melt = melt(FZ.QC.meta.merge.Bottom)
FZ.QC.meta.merge.Bottom.melt$Cluster = gsub("Cluster_", "", FZ.QC.meta.merge.Bottom.melt$PAM)
FZ.QC.meta.merge.Bottom.melt$Cluster = as.factor(FZ.QC.meta.merge.Bottom.melt$Cluster)
FZ.QC.meta.merge.Bottom.melt.Figure = FZ.QC.meta.merge.Bottom.melt
FZ.QC.meta.merge.Bottom.melt.Figure = FZ.QC.meta.merge.Bottom.melt.Figure[FZ.QC.meta.merge.Bottom.melt.Figure$variable %in% c("NonOverFully.Ratio","SemiOverFully.Ratio","Bottom.Ratio"),]
FZ.QC.meta.merge.Bottom.melt.Figure$variable = gsub("NonOverFully.Ratio","NonKR peptides",FZ.QC.meta.merge.Bottom.melt.Figure$variable)
FZ.QC.meta.merge.Bottom.melt.Figure$variable = gsub("SemiOverFully.Ratio","SemiKR peptides",FZ.QC.meta.merge.Bottom.melt.Figure$variable)
FZ.QC.meta.merge.Bottom.melt.Figure$variable = gsub("Bottom.Ratio","Elevated NonKR peptides",FZ.QC.meta.merge.Bottom.melt.Figure$variable)
FZ.QC.meta.merge.Bottom.melt.Figure$variable = factor(FZ.QC.meta.merge.Bottom.melt.Figure$variable, levels = c("SemiKR peptides","NonKR peptides","Elevated NonKR peptides"))

FZ.QC.meta.merge.Bottom.melt.Figure.nonKR = FZ.QC.meta.merge.Bottom.melt.Figure[FZ.QC.meta.merge.Bottom.melt.Figure$variable == "NonKR peptides",]
FZ.QC.meta.merge.Bottom.melt.Figure.SemiKR = FZ.QC.meta.merge.Bottom.melt.Figure[FZ.QC.meta.merge.Bottom.melt.Figure$variable == "SemiKR peptides",]
FZ.QC.meta.merge.Bottom.melt.Figure.MMPnonKR = FZ.QC.meta.merge.Bottom.melt.Figure[FZ.QC.meta.merge.Bottom.melt.Figure$variable == "Elevated NonKR peptides",]

write.table(FZ.QC.meta.merge.Bottom.melt.Figure.nonKR, "FZ.QC.meta.merge.Bottom.melt.Figure.nonKR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(FZ.QC.meta.merge.Bottom.melt.Figure.SemiKR, "FZ.QC.meta.merge.Bottom.melt.Figure.SemiKR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(FZ.QC.meta.merge.Bottom.melt.Figure.MMPnonKR, "FZ.QC.meta.merge.Bottom.melt.Figure.MMPnonKR.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Figure 7D
############################################################
theme_set(theme_grey(base_size = 20))
pdf("FZ_Boxplot_Ratios_Cluster.pdf", width = 6, height = 5)
ggplot(FZ.QC.meta.merge.Bottom.melt.Figure.nonKR, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Cluster), size = 2, shape = 21, position = position_jitterdodge()) +
  theme_bw() + ylab("NonKR ratio") + xlab("clusters") + stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")

ggplot(FZ.QC.meta.merge.Bottom.melt.Figure.SemiKR, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Cluster), size = 2, shape = 21, position = position_jitterdodge()) +
  theme_bw() + ylab("SemiKR ratio") + xlab("clusters") + stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")

ggplot(FZ.QC.meta.merge.Bottom.melt.Figure.MMPnonKR, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Cluster), size = 2, shape = 21, position = position_jitterdodge()) +
  theme_bw() + ylab("NonKR peptides enriched \nin ECM-related proteins") + xlab("clusters") + stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")
dev.off()

############################################################
## Figure S7D
############################################################
FD.meta = read.table("input_metadata/FD_GLBLprot_metadata.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)
FD.protein = read.table("input_protein/FD_GLBLprot_MI_FDbridge_Abund_20201002_Imput_pre_Ave.tsv", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = 1)
FD.cluster.order.1 = c("FD_p10_FH_25_222844; 222854","FD_p11_FH_36_265490; 265491","FD_p20_FH_26_223154; 222912","FD_p01_FH_31_240716","FD_p02_FH_76_493391","FD_p01_FH_57_372440","FD_p16_FH_17_212713","FD_p11_FH_39_271005","FD_p14_FH_53_350407","FD_p15_FH_29_232982","FD_p15_FH_4_206053","FD_p05_FH_84_533557","FD_p03_FH_34_260388","FD_p07_FH_5_500475","FD_p07_FH_83_517774","FD_p06_FH_85_555178; 555177","FD_p18_FH_1_500263","FD_p14_UA_CU-29_S-07-0001549_A3","FD_p08_FH_48_297929","FD_p12_FH_3_500395","FD_p02_FH_11_208877","FD_p02_MC_4086_39149_C6","FD_p03_FH_65_465021; 465023","FD_p17_FH_70_494830; 494829","FD_p02_FH_28_231865; 230368","FD_p08_FH_86_535474","FD_p12_FH_97_352244","FD_p21_UA_CU-18_S-14-0036688_B4","FD_p04_FH_19_215376","FD_p13_FH_35_262427","FD_p05_FH_2_205927","FD_p08_UA_CU-21_S-17-0012504_C1","FD_p07_FH_15_210039","FD_p13_UA_CU-34_S-13-0008803_J1","FD_p13_UA_CU-38_S-13-0021879_A8","FD_p16_UA_CU-24_S-07-0006305_A1","FD_p17_UA_CU-28_S-07-0021920_B3","FD_p20_UA_CU-32_S-13-0010181_B3","FD_p11_UA_CU-46_S-12-0007994_B6","FD_p06_UA_CU-31_S-07-0032072_A8","FD_p06_FH_94_654490","FD_p10_MC_5161_41634_A1","FD_p07_UA_CU-39_S-11-0009142_B2","FD_p19_UA_CU-44_S-12-0006563_D3","FD_p12_UA_CU-40_S-14-0005518_B2","FD_p12_UA_CU-26_S-09-0003181_A5","FD_p04_UA_CU-30_S-10-0035607_A4","FD_p03_UA_CU-35_S-14-0030753_E1","FD_p01_UA_CU-42_S-11-0002981_E3","FD_p10_MC_7924_41550_A1","FD_p16_MC_3219_41837_C2","FD_p21_UA_CU-33_S-14-0021444_B4","FD_p13_UA_CU-27_S-09-0000175_C3","FD_p19_UA_CU-14_S-13-0003847_B3","FD_p01_UA_CU-23_S-15-0033313_B1","FD_p09_MC_2054_36978_C2","FD_p15_UA_CU-02_S-13-0021202_A5","FD_p17_UA_CU-16_S-15-0016317_B7","FD_p06_UA_CU-01_S-07-0009498_A10","FD_p21_MC_3480_38147_A2","FD_p05_UA_CU-17_S-13-0008240_A1","FD_p04_UA_CU-15_S-13-0034287_A2","FD_p20_MC_7339_40414_A1","FD_p07_UA_CU-19_S-17-0024563_A14","FD_p17_MC_4421_39549_D1","FD_p06_MC_1711_36173_D1","FD_p05_MC_4378_39433_A4","FD_p01_MC_7053_39955_B1","FD_p03_FH_95_631766","FD_p07_MC_7715_41079_A6","FD_p19_MC_4114_39184_A1","FD_p04_MC_2313_37530_B9","FD_p21_MC_7203_40170_B2","FD_p15_MC_7078_39996_A2","FD_p09_FH_72_476540","FD_p05_MC_2398_37572_A1","FD_p18_MC_4350_38999_B2","FD_p12_MC_8318_41218_C3","FD_p02_MC_4787_37440_A1")
FD.cluster.order.2 = c("FD_p06_MC_8986_36754_C4","FD_p03_MC_2710_37852_F3","FD_p02_MC_7598_40885_B1","FD_p15_MC_8985_42454_F2","FD_p08_MC_13280_42695_C3","FD_p11_MC_4237_39352_B1","FD_p08_UA_CU-22_S-15-0030998_A1","FD_p11_MC_8005_39918_C2","FD_p16_MC_3812_38735_A6","FD_p21_MC_4322_38538_C14","FD_p17_MC_4968_37624_B1","FD_p14_MC_2303_37504_B1","FD_p03_FH_56_360612","FD_p06_MC_3845_38791_A1","FD_p01_MC_7301_40353_A3","FD_p19_UA_CU-13_S-13-0006146_C3","FD_p14_MC_7064_39968_D7","FD_p15_UA_CU-07_S-08-0003179_A1","FD_p05_MC_7270_40298_K1","FD_p04_FH_96_206843","FD_p01_MC_13136_42571_C1","FD_p17_FH_7_206382; 206375","FD_p19_FH_99_418442; 418441","FD_p04_UA_CU-45_S-11-0025953_B7","FD_p08_UA_CU-36_S-12-0032630_B14","FD_p18_MC_4541_41911_C6","FD_p18_MC_4050_39099_A2","FD_p15_FH_38_265871","FD_p10_UA_CU-43_S-14-0024532_A2","FD_p18_MC_2807_37917_D6","FD_p14_FH_49_296927","FD_p17_UA_CU-41_S-11-0016358_A3","FD_p08_FH_90_555596","FD_p04_UA_CU-25_S-11-0010626_B3","FD_p21_FH_20_218201","FD_p05_MC_3714_38616_A1","FD_p12_UA_CU-20_S-15-0035461_B3","FD_p11_UA_CU-04_S-10-0000363_C1","FD_p02_UA_CU-12_S-14-0007480_A6","FD_p20_MC_2685_37845_B1","FD_p13_MC_13129_42569_E1","FD_p14_UA_CU-03_S-08-0007631_A2","FD_p07_MC_5462_38435_B1","FD_p09_UA_CU-06_S-07-0023499_A2","FD_p07_UA_CU-11_S-13-0022284_A4","FD_p18_UA_CU-10_S-14-0022078_C1","FD_p11_UA_CU-09_S-11-0017293_B1","FD_p09_UA_CU-08_S-07-0027436_B6","FD_p10_MC_2510_37571_A1","FD_p04_FH_22_218317; 218726","FD_p18_FH_69_536848; 533117","FD_p09_FH_41_277568","FD_p16_FH_82_517758","FD_p19_FH_71_474534; 474536","FD_p21_FH_27_228608","FD_p09_FH_78_527499","FD_p16_FH_12_501301","FD_p19_FH_79_503836; 503834","FD_p20_FH_55_351913; 352452","FD_p12_MC_2553_37712_C1","FD_p10_FH_6_500550","FD_p11_FH_77_487495","FD_p13_FH_93_555588","FD_p09_FH_81_525976; 525977","FD_p14_UA_CU-37_S-11-0020094_A3","FD_p01_MC_8228_41080_A1","FD_p08_MC_3604_38471_A2","FD_p20_FH_91_555175","FD_p20_FH_80_476423","FD_p03_FH_33_244726; 246014","FD_p03_FH_14_210413","FD_p10_FH_8_501213","FD_p13_MC_2121_36517_G2","FD_p02_FH_52_319488","FD_p16_FH_10_500869","FD_p06_FH_44_277998","FD_p15_UA_CU-05_S-07-0019837_B1","FD_p13_MC_5980_39738_D2","FD_p05_MC_4564_41946_A1")
FD.cluster.order = c(FD.cluster.order.1, FD.cluster.order.2)

FD.meta.order = FD.meta
rownames(FD.meta.order) = FD.meta.order$Sample
FD.meta.order = FD.meta.order[FD.cluster.order,]

FD.protein.all.order = FD.protein[, FD.cluster.order]
FD.protein.all.order = as.data.frame(t(FD.protein.all.order))
FD.protein.all.order$subtype = FD.meta.order$Cluster
FD.protein.all.order$response = FD.meta.order$Response
FD.protein.all.order$subtype = gsub("Cluster_","",FD.protein.all.order$subtype)
FD.protein.all.order$subtype = factor(FD.protein.all.order$subtype, levels = c(1,2,3,4,5))
FD.protein.all.order$Sample = rownames(FD.protein.all.order)
FD.protein.all.order.melt = melt(FD.protein.all.order)

MMP.path = c("MMP2","MMP7","MMP8","MMP9","MMP11","MMP14","MMP15","MMP24OS","MMP25","TGFB1","TGFB1I1")

FD.MMP.plot = FD.protein.all.order.melt[FD.protein.all.order.melt$variable %in% MMP.path,]
FD.MMP.plot$variable = factor(FD.MMP.plot$variable, levels = MMP.path)

FD.MMP.plot.1 = FD.MMP.plot[FD.MMP.plot$variable == "MMP2",]
FD.MMP.plot.1 <- droplevels(FD.MMP.plot.1)

theme_set(theme_grey(base_size = 20))
pdf("FD_MMP_Plot.pdf", width = 20, height = 8)
ggplot(FD.MMP.plot, aes(x = variable, y = value, fill = subtype)) + geom_boxplot(outlier.shape = NA) +
  theme_bw() + ylab("relative protein abundance") + xlab("clusters") + stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")
dev.off()

############################################################
## Figure S7G
############################################################
theme_set(theme_grey(base_size = 20))
pdf("FD_Boxplot_Ratios_Cluster_Response.pdf", width = 10, height = 5)
ggplot(FD.QC.meta.merge.Bottom.melt.Figure.MMPnonKR, aes(x = Cluster, y = value, fill = Response)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(fill = Response), size = 2, shape = 21, position = position_jitterdodge()) +
  theme_bw() + ylab("NonKR peptides enriched \nin ECM-related proteins") + xlab("clusters") + stat_compare_means(method = "wilcox.test") +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right")
dev.off()
