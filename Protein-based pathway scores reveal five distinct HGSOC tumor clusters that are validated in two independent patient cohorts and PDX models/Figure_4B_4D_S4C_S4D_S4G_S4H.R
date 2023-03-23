############################################################
## title: "Figure 4B_4D_S4C_S4D_S4G_S4H"
## author: "Oscar Murillo"
## date: "March 22, 2023"
## panels: "4B, 4D, S4C, S4D, S4G, S4H"
############################################################
############################################################
## Read in Pathway scores for FFPE Discovery
############################################################
scaleFUN = function(x) sprintf("%.1f", x)

FD.pathway = read.table("input_pathway/FD_pathway_20210927.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)

############################################################
## Estimate cluster centroids
############################################################
FD.pathway.matrix = as.matrix(FD.pathway)
FD.pathway.y = FD.meta$Cluster
FD.pathway.genenames = rownames(FD.pathway)
FD.pathway.geneid = rownames(FD.pathway)
FD.pathway.samples = colnames(FD.pathway)

FD.pathway.input = list(x = FD.pathway.matrix, y = FD.pathway.y, genenames = FD.pathway.genenames, geneid = FD.pathway.geneid, sampleid = FD.pathway.samples)
FD.pathway.training = pamr.train(FD.pathway.input)
FD.pathway.results = pamr.cv(FD.pathway.training, FD.pathway.input)
FD.pathway.threshold = 0

FD.pathway.top = pamr.listgenes(FD.pathway.training, FD.pathway.input, threshold = FD.pathway.threshold)
FD.pathway.top = as.data.frame(FD.pathway.top)
rownames(FD.pathway.top) = FD.pathway.top$id
FD.pathway.top$id = NULL
colnames(FD.pathway.top) = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")
FD.pathway.top$Pathway = rownames(FD.pathway.top)

FD.pathway.top = FD.pathway.top[rev(pathway.order),]

############################################################
## Prepare pathway scores, RNA:protein correlation, metadata
############################################################
pathway.order = rownames(FD.pathway)
FD.cluster.order.1 = c("FD_p10_FH_25_222844; 222854","FD_p11_FH_36_265490; 265491","FD_p20_FH_26_223154; 222912","FD_p01_FH_31_240716","FD_p02_FH_76_493391","FD_p01_FH_57_372440","FD_p16_FH_17_212713","FD_p11_FH_39_271005","FD_p14_FH_53_350407","FD_p15_FH_29_232982","FD_p15_FH_4_206053","FD_p05_FH_84_533557","FD_p03_FH_34_260388","FD_p07_FH_5_500475","FD_p07_FH_83_517774","FD_p06_FH_85_555178; 555177","FD_p18_FH_1_500263","FD_p14_UA_CU-29_S-07-0001549_A3","FD_p08_FH_48_297929","FD_p12_FH_3_500395","FD_p02_FH_11_208877","FD_p02_MC_4086_39149_C6","FD_p03_FH_65_465021; 465023","FD_p17_FH_70_494830; 494829","FD_p02_FH_28_231865; 230368","FD_p08_FH_86_535474","FD_p12_FH_97_352244","FD_p21_UA_CU-18_S-14-0036688_B4","FD_p04_FH_19_215376","FD_p13_FH_35_262427","FD_p05_FH_2_205927","FD_p08_UA_CU-21_S-17-0012504_C1","FD_p07_FH_15_210039","FD_p13_UA_CU-34_S-13-0008803_J1","FD_p13_UA_CU-38_S-13-0021879_A8","FD_p16_UA_CU-24_S-07-0006305_A1","FD_p17_UA_CU-28_S-07-0021920_B3","FD_p20_UA_CU-32_S-13-0010181_B3","FD_p11_UA_CU-46_S-12-0007994_B6","FD_p06_UA_CU-31_S-07-0032072_A8","FD_p06_FH_94_654490","FD_p10_MC_5161_41634_A1","FD_p07_UA_CU-39_S-11-0009142_B2","FD_p19_UA_CU-44_S-12-0006563_D3","FD_p12_UA_CU-40_S-14-0005518_B2","FD_p12_UA_CU-26_S-09-0003181_A5","FD_p04_UA_CU-30_S-10-0035607_A4","FD_p03_UA_CU-35_S-14-0030753_E1","FD_p01_UA_CU-42_S-11-0002981_E3","FD_p10_MC_7924_41550_A1","FD_p16_MC_3219_41837_C2","FD_p21_UA_CU-33_S-14-0021444_B4","FD_p13_UA_CU-27_S-09-0000175_C3","FD_p19_UA_CU-14_S-13-0003847_B3","FD_p01_UA_CU-23_S-15-0033313_B1","FD_p09_MC_2054_36978_C2","FD_p15_UA_CU-02_S-13-0021202_A5","FD_p17_UA_CU-16_S-15-0016317_B7","FD_p06_UA_CU-01_S-07-0009498_A10","FD_p21_MC_3480_38147_A2","FD_p05_UA_CU-17_S-13-0008240_A1","FD_p04_UA_CU-15_S-13-0034287_A2","FD_p20_MC_7339_40414_A1","FD_p07_UA_CU-19_S-17-0024563_A14","FD_p17_MC_4421_39549_D1","FD_p06_MC_1711_36173_D1","FD_p05_MC_4378_39433_A4","FD_p01_MC_7053_39955_B1","FD_p03_FH_95_631766","FD_p07_MC_7715_41079_A6","FD_p19_MC_4114_39184_A1","FD_p04_MC_2313_37530_B9","FD_p21_MC_7203_40170_B2","FD_p15_MC_7078_39996_A2","FD_p09_FH_72_476540","FD_p05_MC_2398_37572_A1","FD_p18_MC_4350_38999_B2","FD_p12_MC_8318_41218_C3","FD_p02_MC_4787_37440_A1")
FD.cluster.order.2 = c("FD_p06_MC_8986_36754_C4","FD_p03_MC_2710_37852_F3","FD_p02_MC_7598_40885_B1","FD_p15_MC_8985_42454_F2","FD_p08_MC_13280_42695_C3","FD_p11_MC_4237_39352_B1","FD_p08_UA_CU-22_S-15-0030998_A1","FD_p11_MC_8005_39918_C2","FD_p16_MC_3812_38735_A6","FD_p21_MC_4322_38538_C14","FD_p17_MC_4968_37624_B1","FD_p14_MC_2303_37504_B1","FD_p03_FH_56_360612","FD_p06_MC_3845_38791_A1","FD_p01_MC_7301_40353_A3","FD_p19_UA_CU-13_S-13-0006146_C3","FD_p14_MC_7064_39968_D7","FD_p15_UA_CU-07_S-08-0003179_A1","FD_p05_MC_7270_40298_K1","FD_p04_FH_96_206843","FD_p01_MC_13136_42571_C1","FD_p17_FH_7_206382; 206375","FD_p19_FH_99_418442; 418441","FD_p04_UA_CU-45_S-11-0025953_B7","FD_p08_UA_CU-36_S-12-0032630_B14","FD_p18_MC_4541_41911_C6","FD_p18_MC_4050_39099_A2","FD_p15_FH_38_265871","FD_p10_UA_CU-43_S-14-0024532_A2","FD_p18_MC_2807_37917_D6","FD_p14_FH_49_296927","FD_p17_UA_CU-41_S-11-0016358_A3","FD_p08_FH_90_555596","FD_p04_UA_CU-25_S-11-0010626_B3","FD_p21_FH_20_218201","FD_p05_MC_3714_38616_A1","FD_p12_UA_CU-20_S-15-0035461_B3","FD_p11_UA_CU-04_S-10-0000363_C1","FD_p02_UA_CU-12_S-14-0007480_A6","FD_p20_MC_2685_37845_B1","FD_p13_MC_13129_42569_E1","FD_p14_UA_CU-03_S-08-0007631_A2","FD_p07_MC_5462_38435_B1","FD_p09_UA_CU-06_S-07-0023499_A2","FD_p07_UA_CU-11_S-13-0022284_A4","FD_p18_UA_CU-10_S-14-0022078_C1","FD_p11_UA_CU-09_S-11-0017293_B1","FD_p09_UA_CU-08_S-07-0027436_B6","FD_p10_MC_2510_37571_A1","FD_p04_FH_22_218317; 218726","FD_p18_FH_69_536848; 533117","FD_p09_FH_41_277568","FD_p16_FH_82_517758","FD_p19_FH_71_474534; 474536","FD_p21_FH_27_228608","FD_p09_FH_78_527499","FD_p16_FH_12_501301","FD_p19_FH_79_503836; 503834","FD_p20_FH_55_351913; 352452","FD_p12_MC_2553_37712_C1","FD_p10_FH_6_500550","FD_p11_FH_77_487495","FD_p13_FH_93_555588","FD_p09_FH_81_525976; 525977","FD_p14_UA_CU-37_S-11-0020094_A3","FD_p01_MC_8228_41080_A1","FD_p08_MC_3604_38471_A2","FD_p20_FH_91_555175","FD_p20_FH_80_476423","FD_p03_FH_33_244726; 246014","FD_p03_FH_14_210413","FD_p10_FH_8_501213","FD_p13_MC_2121_36517_G2","FD_p02_FH_52_319488","FD_p16_FH_10_500869","FD_p06_FH_44_277998","FD_p15_UA_CU-05_S-07-0019837_B1","FD_p13_MC_5980_39738_D2","FD_p05_MC_4564_41946_A1")
FD.cluster.order = c(FD.cluster.order.1, FD.cluster.order.2)

FD.pathway.order = FD.pathway[, FD.cluster.order]
FD.meta.order = FD.meta
rownames(FD.meta.order) = FD.meta.order$Sample
FD.meta.order = FD.meta.order[FD.cluster.order,]

FD.pathway.order.t.scale = scale(t(FD.pathway.order))
FD.pathway.order.t.scale[is.na(FD.pathway.order.t.scale)] <- 0
FD.pathway.order.t.scale.t = as.matrix(t(FD.pathway.order.t.scale))

FD.order.clusters = FD.meta.order$Cluster
FD.cluster.pca = as.data.frame(t(FD.pathway.order))
FD.cluster.label = substr(rownames(FD.cluster.pca), start = 1, stop = 2)
FD.cluster.id = FD.order.clusters

## Generate RNA:Protein correlation bar plot
FD.RNA.pathway = read.table("input_pathway/FD_RNA_pathway_20220208.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)

FD.pathway.order.RNAoverlap = FD.pathway.order[,colnames(FD.RNA.pathway)]
FD.pathway.order.RNAoverlap = as.data.frame(t(FD.pathway.order.RNAoverlap))
FD.RNA.pathway = as.data.frame(t(FD.RNA.pathway))

FD.Pro.RNA.pathway = cor(FD.pathway.order.RNAoverlap, FD.RNA.pathway)
FD.Pro.RNA.pathway.df = as.data.frame(diag(FD.Pro.RNA.pathway))
colnames(FD.Pro.RNA.pathway.df) = "protein_RNA"
FD.Figure3.meta.Row.Right = rowAnnotation("corr" = anno_barplot(FD.Pro.RNA.pathway.df[,1], baseline = 0, gp = gpar(fill = c("darkgrey"), col = c("darkgrey")), border = FALSE))

FD.Association = c(5.38,8.60,5.25,11.51,3.69,11.10,4.61,10.43,9.02,8.32,5.62,12.27,10.82,4.67,4.60,4.76,9.42,7.74,5.38,4.67,4.67,5.44,5.46,12.02,10.19,7.81,4.73,2.26,2.89,2.94,2.34,2.62,4.22,2.16,2.69,4.67,3.33,2.70,3.27,2.15,2.45,6.42,2.07,2.45,2.77,3.33,2.41,2.83,2.60,2.26,2.30,2.64,2.60,2.35,2.62,3.05,13.72,3.52,6.02,6.97,11.50,3.65,3.53,4.88,13.86,11.50,15.34,15.42,6.60,14.66,8.72,11.03,14.98,19.33,9.35,7.43,10.51,13.68,10.51,15.47,8.56,14.66,6.73,11.46,12.27,14.66,8.56,6.57,14.88,19.37,19.39,9.17,8.06,9.87,5.06,2.12,2.22,2.87,2.03,2.12,2.89,7.88,2.24,3.80,4.02,2.24,3.67,2.79,4.02,4.67,2.96,2.43,2.32,2.15,2.24,3.10,3.02,3.97,3.83,6.10,3.49,3.21,11.56,18.56,2.96,8.91,6.46,12.26,7.65,13.72,9.50,12.85,11.10,8.04,2.03,21.23,8.12,5.67,2.87,12.77,2.89,2.24,2.38,2.28,3.63,15.42,19.06,3.18,18.90,6.60)
FD.Figure3.Association.Right = rowAnnotation("FDR" = anno_barplot(FD.Association, baseline = 0, gp = gpar(fill = c("black"), col = c("black")), border = FALSE))

## ComplexHeatmap
FD.pathway.top.HM = FD.pathway.top[pathway.order,]
FD.pathway.top.HM$Pathway = NULL
FD.pathway.top.HM = as.data.frame(FD.pathway.top.HM)

FD.pathway.top.HM = matrix(as.numeric(unlist(FD.pathway.top.HM)), ncol = ncol(FD.pathway.top.HM))
colnames(FD.pathway.top.HM) = c("1","2","3","4","5")
FD.pathway.top.HM = as.data.frame(FD.pathway.top.HM)

FD.Figure3.meta.Row = rowAnnotation("1" = anno_barplot(FD.pathway.top.HM[,1], baseline = 0, gp = gpar(fill = c("#F8766D"), col = c("#F8766D")), border = FALSE),
                                    "2" = anno_barplot(FD.pathway.top.HM[,2], baseline = 0, gp = gpar(fill = c("#BB9D00"), col = c("#BB9D00")), border = FALSE),
                                    "3" = anno_barplot(FD.pathway.top.HM[,3], baseline = 0, gp = gpar(fill = c("#00B81F"), col = c("#00B81F")), border = FALSE),
                                    "4" = anno_barplot(FD.pathway.top.HM[,4], baseline = 0, gp = gpar(fill = c("#00A5FF"), col = c("#00A5FF")), border = FALSE),
                                    "5" = anno_barplot(FD.pathway.top.HM[,5], baseline = 0, gp = gpar(fill = c("#E76BF3"), col = c("#E76BF3")), border = FALSE))

FD.Figure3.meta = read.table("Figure_3/FD_metadata.txt", sep = "\t", header = TRUE)
rownames(FD.Figure3.meta) = FD.Figure3.meta$Sample
FD.Figure3.meta.order = FD.Figure3.meta[FD.cluster.order,]
FD.Figure3.meta.order$subtype = factor(FD.Figure3.meta.order$subtype, levels = c(1,2,3,4,5))
FD.Figure3.meta.order$RNA_cluster = factor(FD.Figure3.meta.order$RNA_cluster, levels = c(1,2,3,4,5))
FD.Figure3.meta.order$location = ifelse(FD.Figure3.meta.order$Sample == "FD_p05_MC_4564_41946_A1", "other", FD.Figure3.meta.order$location)
FD.Figure3.meta.order$location = ifelse(FD.Figure3.meta.order$Sample == "FD_p20_MC_2685_37845_B1", "other", FD.Figure3.meta.order$location)

FD.sample.age.col = colorRamp2(c(0, 10, 20), c("white", "lightblue", "blue"))
FD.tumor.purity.col = colorRamp2(c(0.6, 0.8, 1), c("white", "orange", "red"))

FD.protein.VerhakkSubtypes.table$Sample = rownames(FD.protein.VerhakkSubtypes.table)

FD.Figure3.meta.order.1 = merge(FD.Figure3.meta.order, FD.protein.VerhakkSubtypes.table, by = "Sample")
FD.Figure3.meta.order.1$Subtype = gsub("PRO","Proliferative", FD.Figure3.meta.order.1$Subtype)
FD.Figure3.meta.order.1$Subtype = gsub("MES","Mesenchymal", FD.Figure3.meta.order.1$Subtype)
FD.Figure3.meta.order.1$Subtype = gsub("IMR","Immunoreactive", FD.Figure3.meta.order.1$Subtype)
FD.Figure3.meta.order.1$Subtype = gsub("DIF","Differentiated", FD.Figure3.meta.order.1$Subtype)
rownames(FD.Figure3.meta.order.1) = FD.Figure3.meta.order.1$Sample
FD.Figure3.meta.order.1 = FD.Figure3.meta.order.1[FD.Figure3.meta.order$Sample,]

FD.Figure3.meta.Ann = HeatmapAnnotation(site = FD.Figure3.meta.order.1$site, sample.age = FD.Figure3.meta.order.1$age, location = FD.Figure3.meta.order.1$location,
                                        neo.adjuvant = FD.Figure3.meta.order.1$neo.adjuvant, TCGA.subtype = FD.Figure3.meta.order.1$Subtype, tumor.purity = FD.Figure3.meta.order.1$tumor.purity,
                                        genome.instability = FD.Figure3.meta.order.1$genome.instability,
                                        response = FD.Figure3.meta.order.1$response, subtype = FD.Figure3.meta.order.1$subtype,
                                        col = list(site = c("FHCRC" = "turquoise1","Mayo" = "blue","UAB" = "green"),
                                                   location = c("OV" = "yellow","OM" = "red","other" = "dodgerblue","mix" = "orange"),
                                                   response = c("sensitive" = "gold","refractory" = "blue"),
                                                   neo.adjuvant = c("yes" = "black","no" = "white"),
                                                   sample.age = FD.sample.age.col, tumor.purity = FD.tumor.purity.col,
                                                   TCGA.subtype = c("Proliferative" = "blue","Mesenchymal" = "violet","Immunoreactive" = "red","Differentiated" = "green"),
                                                   subtype = c("1" = "#F8766D","2" = "#BB9D00","3" = "#00B81F","4" = "#00A5FF","5" = "#E76BF3")),
                                        na_col = "lightgray")

FD.Figure3.pathwayGroups = read.table("Figure_3/FD_Pathway_Groups.txt", sep = "\t", header = TRUE)
FD.label.pathways.2 = c("","","","","","","","","","","Complement cascade","","","","","Role of phospholipids in phagocytosis","","","","","","","","","",
                        "","","","Antigen processing and presentation","","","Interferon alpha response","","","Interferon gamma response","","","Adaptive immune system","","","","Coagulation","","","Epithelial mesenchymal transition","","","Heme metabolism","","","Hypoxia","","","TGFb signature","","","","","","Ribosome","","","","Developmental biology","","","","Eukaryotic translation elongation","","","","Eukaryotic translation initiation","","","","Nonsense mediated decay (NMD)","","","","Response of EIF2AK4 (GCN2) \nto amino acid deficiency","","","","Selenoamino acid metabolism","","","","rRNA processing","","","","","DNA replication","","","","E2F targets","","","","G2M checkpoint","","","","DNA replication","","","","DNA strand elongation","","","","Metabolism of RNA","","","","Unwinding of DNA","","","","Core mitotic","","","","","","Adipogenesis","","","","","Oxidative phosphorylation","","","","","Metabolism","","","","","Complex I biogenesis","","","","","The citric acid(TCA) cycle and\n respiratory electron transport","","","")

############################################################
## Figure 4B
############################################################
theme_set(theme_grey(base_size = 16))
pdf("FD_Pathway_Heatmap.pdf", height = 14, width = 16, useDingbats = FALSE)
ht_list = Heatmap(FD.pathway.order.t.scale.t, top_annotation = FD.Figure3.meta.Ann,
                  cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(FD.Figure3.meta.order.1$subtype, levels = c("1","2","3","4","5")),
                  row_split = factor(FD.Figure3.pathwayGroups$Group, levels = c("6","5","4","3","2","1")), border = TRUE, 
                  row_title = c("","","","","",""), column_title = c("","","","",""),
                  show_column_names = FALSE, row_labels = FD.label.pathways.2, left_annotation = FD.Figure3.meta.Row, right_annotation = FD.Figure3.Association.Right)
draw(ht_list, heatmap_legend_side = c("bottom"), annotation_legend_side = "bottom",  padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

############################################################
## Figure 4D
############################################################
FD.LOH.data = read.table("LOH.txt", header = TRUE, sep = "\t", quote = NULL, row.names = NULL)
FD.LOH.data$chr17LOH = factor(FD.LOH.data$chr17LOH)
FD.LOH.data = FD.LOH.data[,c("Sample","chr17LOH")]
FD.Figure3.meta.order.2 = merge(FD.Figure3.meta.order.1, FD.LOH.data, by = "Sample", all.x = TRUE)
FD.Figure3.meta.order.2.LOH = FD.Figure3.meta.order.2[!is.na(FD.Figure3.meta.order.2$chr17LOH),]

LOH.bar = ggplot(data = FD.Figure3.meta.order.2.LOH, aes(x = subtype, y = ..count../ sum(..count..), fill = chr17LOH)) + 
  geom_bar(position = "fill") + scale_fill_manual("legend", values = c("1" = "black", "0" = "white"), na.value = "gray") +
  xlab("cluster") + ylab("% of samples") + theme_bw() + scale_y_continuous(labels = scales::percent_format()) + theme(legend.position = "none") + 
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank())

############################################################
## Figure S4C
############################################################
response.bar = ggplot(data = FD.Figure3.meta.order.1, aes(x = subtype, y = ..count../ sum(..count..), fill = response)) + 
  geom_bar(position = "fill") + scale_fill_manual("legend", values = c("sensitive" = "gold","refractory" = "blue"), na.value = "gray") +
  xlab("cluster") + ylab("% of samples") + theme_bw() + scale_y_continuous(labels = scales::percent_format()) + theme(legend.position = "none") + 
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank())

############################################################
## Figure S4D
############################################################
location.bar = ggplot(data = FD.Figure3.meta.order.1, aes(x = subtype, y = ..count../ sum(..count..), fill = location)) + 
  geom_bar(position = "fill") + scale_fill_manual("legend", values = c("OV" = "yellow","OM" = "red","other" = "dodgerblue","mix" = "orange"), na.value = "gray") +
  ylab("tumor location") + theme_bw() + scale_y_continuous(labels = scales::percent_format()) + theme(legend.position = "none") + 
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank())

TCGA.bar = ggplot(data = FD.Figure3.meta.order.1, aes(x = subtype, y = ..count../ sum(..count..), fill = Subtype)) + 
  geom_bar(position = "fill") + scale_fill_manual("legend", values = c("Proliferative" = "blue","Mesenchymal" = "violet","Immunoreactive" = "red","Differentiated" = "green"), na.value = "gray") +
  ylab("TCGA.subtype") + theme_bw() + scale_y_continuous(labels = scales::percent_format()) + theme(legend.position = "none") + 
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank())

############################################################
## Set of analysis for Figure S4G
############################################################
library(org.Hs.eg.db)
hs = org.Hs.eg.db

OV.Beta.Sub = OV.Beta
OV.Beta.Sub$SYMBOL = rownames(OV.Beta.Sub)

newGeneId = select(hs, keys = rownames(OV.Beta.Sub), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
OV.Beta.Sub.NewGeneID = merge(OV.Beta.Sub, newGeneId, by = "SYMBOL")
OV.Beta.Sub.NewGeneID = na.omit(OV.Beta.Sub.NewGeneID)
rownames(OV.Beta.Sub.NewGeneID) = OV.Beta.Sub.NewGeneID$ENTREZID
OV.Beta.Sub.NewGeneID$ENTREZID = NULL
OV.Beta.Sub.NewGeneID$SYMBOL = NULL

OV.Beta.VerhakkSubtypes = get.subtypes(as.matrix(OV.Beta.Sub.NewGeneID), rownames(OV.Beta.Sub.NewGeneID), method = "Verhaak")
OV.Beta.VerhaakList = as.vector(OV.Beta.VerhakkSubtypes$Verhaak.subtypes)

OV.Beta.VerhakkSubtypes.table = as.data.frame(OV.Beta.VerhakkSubtypes$gsva.out)
OV.Beta.VerhakkSubtypes.table$Subtype = OV.Beta.VerhaakList

write.table(OV.Beta.VerhakkSubtypes.table, file = "FD_RNA_Subtype_20230302.txt", sep = "\t", quote = FALSE, row.names = TRUE)

FD.protein.Sub = FD.protein
FD.protein.Sub$SYMBOL = rownames(FD.protein.Sub)

newGeneId = select(hs, keys = rownames(FD.protein.Sub), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
FD.protein.Sub.NewGeneID = merge(FD.protein.Sub, newGeneId, by = "SYMBOL")
FD.protein.Sub.NewGeneID = na.omit(FD.protein.Sub.NewGeneID)
rownames(FD.protein.Sub.NewGeneID) = FD.protein.Sub.NewGeneID$ENTREZID
FD.protein.Sub.NewGeneID$ENTREZID = NULL
FD.protein.Sub.NewGeneID$SYMBOL = NULL

FD.protein.VerhakkSubtypes = get.subtypes(as.matrix(FD.protein.Sub.NewGeneID), rownames(FD.protein.Sub.NewGeneID), method = "Verhaak")
FD.protein.VerhaakList = as.vector(FD.protein.VerhakkSubtypes$Verhaak.subtypes)

FD.protein.VerhakkSubtypes.table = as.data.frame(FD.protein.VerhakkSubtypes$gsva.out)
FD.protein.VerhakkSubtypes.table$Subtype = FD.protein.VerhaakList
FD.protein.VerhakkSubtypes.table$Sample = rownames(FD.protein.VerhakkSubtypes.table)
FD.protein.VerhakkSubtypes.table = FD.protein.VerhakkSubtypes.table[rownames(FD.Figure3.meta),]
FD.protein.VerhakkSubtypes.table$Sample = NULL

write.table(FD.protein.VerhakkSubtypes.table, file = "FD_Protein_Subtype_20230302.txt", sep = "\t", quote = FALSE, row.names = TRUE)

## Differential Analysis FD RNA Samples based on TCGA Subtypes
tTestOneVsRest.TCGA = perform_t_tests_all_classes_one_vs_rest(dataMatrix = OV.Beta, classVector = OV.Beta.VerhaakList)

TCGA.DiffMeans = as.data.frame(tTestOneVsRest.TCGA$Difference.Between.Means)
colnames(TCGA.DiffMeans) = paste(colnames(TCGA.DiffMeans), "DiffMeans", sep = "_")
TCGA.pValue = as.data.frame(tTestOneVsRest.TCGA$P.Values)
colnames(TCGA.pValue) = paste(colnames(TCGA.pValue), "pValue", sep = "_")

OV.RNA.Diff.Results = cbind(TCGA.DiffMeans, TCGA.pValue)
write.table(OV.RNA.Diff.Results, file = "Pathway/FD_RNA_Diff_Results_20220225.txt", sep = "\t", quote = FALSE, row.names = TRUE)

OV.RNA.Diff.Results$PRO_FDR = p.adjust(OV.RNA.Diff.Results$PRO_pValue, method = "fdr")
OV.RNA.Diff.Results$MES_FDR = p.adjust(OV.RNA.Diff.Results$MES_pValue, method = "fdr")
OV.RNA.Diff.Results$DIF_FDR = p.adjust(OV.RNA.Diff.Results$DIF_pValue, method = "fdr")
OV.RNA.Diff.Results$IMR_FDR = p.adjust(OV.RNA.Diff.Results$IMR_pValue, method = "fdr")
OV.RNA.Diff.Results$Gene = rownames(OV.RNA.Diff.Results)

## Differential Analysis FD RNA Samples based on TCGA Subtypes - Heatmaps
chosenProbes.one.TCGA = c()
for(i in 1:4){
  sigProbes.TCGA = row.names(tTestOneVsRest.TCGA$P.Values[tTestOneVsRest.TCGA$P.Values[,i] < 0.025,])
  sigDiff.TCGA = tTestOneVsRest.TCGA$Difference.Between.Means[which(row.names(tTestOneVsRest.TCGA$P.Values) %in% sigProbes.TCGA), i]
  highDiffProbes.TCGA = names(tail(sort(sigDiff.TCGA), 100))
  print(length(highDiffProbes.TCGA))
  chosenProbes.one.TCGA = c(chosenProbes.one.TCGA, highDiffProbes.TCGA)
}

OV.Beta.t.scale = scale(t(OV.Beta[chosenProbes.one.TCGA,]))
OV.Beta.t.scale[is.na(OV.Beta.t.scale)] <- 0
OV.Beta.t.scale.t = as.matrix(t(OV.Beta.t.scale))

FD.RNA.Ann = HeatmapAnnotation(TCGA.subtype = OV.Beta.VerhaakList,
                               col = list(TCGA.subtype = c("PRO" = "blue","MES" = "violet","IMR" = "red","DIF" = "green")), na_col = "lightgray")

theme_set(theme_grey(base_size = 16))
pdf("Pathway/FD_RNA_TCGA_DiffExp_Heatmap_Top400.pdf", height = 10, width = 12, useDingbats = FALSE)
ht_list = Heatmap(OV.Beta.t.scale.t, top_annotation = FD.RNA.Ann, use_raster = TRUE,
                  cluster_rows = TRUE, cluster_columns = TRUE, column_split = factor(OV.Beta.meta$Verhaak),
                  border = TRUE, show_column_names = FALSE, show_row_names = FALSE)
draw(ht_list, heatmap_legend_side = c("right"), annotation_legend_side = "right",  padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

OV.RNA.Diff.pValue.Cutoff = 0.01
OV.RNA.Diff.Results.PRO = OV.RNA.Diff.Results[OV.RNA.Diff.Results$PRO_pValue <= OV.RNA.Diff.pValue.Cutoff,]
OV.RNA.Diff.Results.MES = OV.RNA.Diff.Results[OV.RNA.Diff.Results$MES_pValue <= OV.RNA.Diff.pValue.Cutoff,]
OV.RNA.Diff.Results.DIF = OV.RNA.Diff.Results[OV.RNA.Diff.Results$DIF_pValue <= OV.RNA.Diff.pValue.Cutoff,]
OV.RNA.Diff.Results.IMR = OV.RNA.Diff.Results[OV.RNA.Diff.Results$IMR_pValue <= OV.RNA.Diff.pValue.Cutoff,]

## PRO = 2318
length(OV.RNA.Diff.Results.PRO$Gene)
## MES = 2808
length(OV.RNA.Diff.Results.MES$Gene)
## DIF = 2606
length(OV.RNA.Diff.Results.DIF$Gene)
## IMR = 2463
length(OV.RNA.Diff.Results.IMR$Gene)

## All 7257
All.Significant = unique(c(OV.RNA.Diff.Results.PRO$Gene, OV.RNA.Diff.Results.MES$Gene, OV.RNA.Diff.Results.DIF$Gene, OV.RNA.Diff.Results.IMR$Gene))
OV.Beta.Significant = OV.Beta[All.Significant,]

OV.Beta.t.scale.2 = scale(t(OV.Beta[All.Significant,]))
OV.Beta.t.scale.2[is.na(OV.Beta.t.scale.2)] <- 0
OV.Beta.t.scale.2.t = as.matrix(t(OV.Beta.t.scale.2))
OV.Beta.t.scale.2.t = as.data.frame(OV.Beta.t.scale.2.t)

OV.Beta.t.scale.2.t[] <- lapply(OV.Beta.t.scale.2.t, function(x) ifelse(x > 4, 4, x))
OV.Beta.t.scale.2.t[] <- lapply(OV.Beta.t.scale.2.t, function(x) ifelse(x < -4, -4, x))
OV.Beta.t.scale.2.t = as.matrix(OV.Beta.t.scale.2.t)

FD.RNA.Ann.2 = HeatmapAnnotation(show_legend = c(FALSE), TCGA.subtype = OV.Beta.VerhaakList,
                                 col = list(TCGA.subtype = c("PRO" = "blue","MES" = "violet","IMR" = "red","DIF" = "green")), na_col = "lightgray")

Nature.Genes = c("FAP","ANGPTL1","MCM2","CXCL10","CXCL11","SOX11","MUC1","HMGA2","ANGPTL2","CXCR3","MUC16","ASLP1","PCNA")
Nature.Genes = intersect(Nature.Genes, rownames(OV.Beta.t.scale.2.t))

Genes.All.reOrg = intersect(rownames(OV.Beta.t.scale.2.t), Nature.Genes)
Genes.All.inx = which(rownames(OV.Beta.t.scale.2.t) %in% Genes.All.reOrg)

FD.RNA.All = rowAnnotation(foo = anno_mark(at = Genes.All.inx, labels = Genes.All.reOrg, side = "left"))

FD.genes.matrix = as.matrix(OV.Beta[All.Significant,])

FD.genes.input = list(x = FD.genes.matrix, y = OV.Beta.VerhaakList, genenames = rownames(FD.genes.matrix), geneid = rownames(FD.genes.matrix), sampleid = colnames(FD.genes.matrix))
FD.genes.training = pamr.train(FD.genes.input)

FD.genes.threshold = 0

FD.genes.top = pamr.listgenes(FD.genes.training, FD.genes.input, threshold = FD.genes.threshold)
FD.genes.top = as.data.frame(FD.genes.top)
rownames(FD.genes.top) = FD.genes.top$id
FD.genes.top$id = NULL
colnames(FD.genes.top) = c("DIF","IMR","MES","PRO")
FD.genes.top$Genes = rownames(FD.genes.top)

FD.genes.top.HM = FD.genes.top
FD.genes.top.HM$Genes = NULL
FD.genes.top.HM = as.data.frame(FD.genes.top.HM[rownames(OV.Beta.t.scale.2.t),])

FD.genes.top.HM.Nature = FD.genes.top.HM[rownames(FD.genes.top.HM) %in% Nature.Genes,]
FD.genes.top.HM.NonNature = FD.genes.top.HM[!rownames(FD.genes.top.HM) %in% Nature.Genes,]
FD.genes.top.HM.NonNature$DIF = 0
FD.genes.top.HM.NonNature$IMR = 0
FD.genes.top.HM.NonNature$MES = 0
FD.genes.top.HM.NonNature$PRO = 0
FD.genes.top.HM = rbind(FD.genes.top.HM.Nature, FD.genes.top.HM.NonNature)
FD.genes.top.HM = as.data.frame(FD.genes.top.HM[rownames(OV.Beta.t.scale.2.t),])

FD.genes.top.HM = matrix(as.numeric(unlist(FD.genes.top.HM)), ncol = ncol(FD.genes.top.HM))
colnames(FD.genes.top.HM) = c("DIF","IMR","MES","PRO")
FD.genes.top.HM = as.data.frame(FD.genes.top.HM)

FD.Genes.RowAnnotation = rowAnnotation("MES" = anno_barplot(FD.genes.top.HM[,3], baseline = 0, ylim = c(-0.6,0.6), gp = gpar(fill = c("violet"), col = c("violet")), border = FALSE, axis_param = list(gp = gpar(fontsize = 10))),
                                       "IMR" = anno_barplot(FD.genes.top.HM[,2], baseline = 0, ylim = c(-0.6,0.6), gp = gpar(fill = c("red"), col = c("red")), border = FALSE, axis_param = list(gp = gpar(fontsize = 10))),
                                       "PRO" = anno_barplot(FD.genes.top.HM[,4], baseline = 0, ylim = c(-0.6,0.6), gp = gpar(fill = c("blue"), col = c("blue")), border = FALSE, axis_param = list(gp = gpar(fontsize = 10))),
                                       "DIF" = anno_barplot(FD.genes.top.HM[,1], baseline = 0, ylim = c(-0.6,0.6), gp = gpar(fill = c("green"), col = c("green")), border = FALSE, axis_param = list(gp = gpar(fontsize = 10))))


text_list = list(text1 = "cytokine signaling",
                 text2 = "ECM interaction",
                 text3 = "DNA replication",
                 text4 = "",
                 text5 = "DNA replication\nmetabolism")
Right_Annotation_Test = rowAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(unlist(text_list)) + unit(4, "mm")))

set.seed(12345)
pdf("Pathway/FD_RNA_TCGA_DiffExp_Heatmap_AllSign_Cluster5.pdf", height = 6, width = 10, useDingbats = FALSE)
ht_list1 = Heatmap(OV.Beta.t.scale.2.t, name = "Z-score", top_annotation = FD.RNA.Ann.2, use_raster = TRUE,
                   cluster_rows = TRUE, cluster_columns = TRUE, column_split = factor(OV.Beta.meta$Verhaak),
                   border = TRUE, show_column_names = FALSE, show_row_names = FALSE, show_row_dend = FALSE,
                   row_names_gp = gpar(fontsize = 32), row_km = 5, row_title = c("","","","",""), 
                   left_annotation = c(FD.RNA.All, FD.Genes.RowAnnotation))
ht_list1 = draw(ht_list1, heatmap_legend_side = c("right"), padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

ht_list2_rows = row_order(ht_list1)

OV.Beta.t.scale.2.t.Group1 = OV.Beta.t.scale.2.t[ht_list2_rows$'1',]
OV.Beta.t.scale.2.t.Group2 = OV.Beta.t.scale.2.t[ht_list2_rows$'2',]
OV.Beta.t.scale.2.t.Group3 = OV.Beta.t.scale.2.t[ht_list2_rows$'5',]
OV.Beta.t.scale.2.t.Group4 = OV.Beta.t.scale.2.t[ht_list2_rows$'3',]
OV.Beta.t.scale.2.t.Group5 = OV.Beta.t.scale.2.t[ht_list2_rows$'4',]

Group1.Genes = rownames(OV.Beta.t.scale.2.t.Group1)
Group2.Genes = rownames(OV.Beta.t.scale.2.t.Group2)
Group3.Genes = rownames(OV.Beta.t.scale.2.t.Group3)
Group4.Genes = rownames(OV.Beta.t.scale.2.t.Group4)
Group5.Genes = rownames(OV.Beta.t.scale.2.t.Group5)

Groups.Names = c("Group1","Group2","Group3","Group4","Group5")
Groups.List = list(Group1.Genes, Group2.Genes, Group3.Genes, Group4.Genes, Group5.Genes)

## ECM
ECM.Paths = c("Extracellular matrix organization","Muscle contraction","Smooth Muscle Contraction","Collagen biosynthesis and modifying enzymes","Axon guidance","Assembly of collagen fibrils and other multimeric structures")

ECM = list()
for(i in ECM.Paths){
  count = 1
  for(j in Groups.List){
    gene.in = Pathways.Reactome[[i]]
    print(gene.in[1:10])
    gene.out = setdiff(rownames(OV.Beta), gene.in)
    
    strata1.in = length(intersect(gene.in, j))
    no_strata1.in = length(setdiff(gene.in, strata1.in))  
    strata1.out = length(intersect(gene.out, j))  
    no_strata1.out = length(setdiff(gene.out, strata1.out))    
    
    fish =  matrix(c(strata1.in, no_strata1.in, strata1.out, no_strata1.out), nrow = 2)    
    test = fisher.test(fish, y = NULL, alternative = "greater")  
    p.val = test$p.value     
    print(p.val)
    ECM[[i]][[count]] = p.val
    count = count + 1
  }
}

ECM.Results = as.data.frame(matrix(unlist(ECM), ncol = length(ECM.Paths)))
colnames(ECM.Results) = ECM.Paths
rownames(ECM.Results) = Groups.Names

## Metabolism
Metabolism.Paths = c("Metabolism","The citric acid (TCA) cycle and respiratory electron transport","Pyruvate metabolism and Citric Acid (TCA) cycle","Citric acid cycle (TCA cycle)","Gluconeogenesis")

Metabolism = list()
for(i in Metabolism.Paths){
  count = 1
  for(j in Groups.List){
    gene.in = Pathways.Reactome[[i]]
    print(gene.in[1:10])
    gene.out = setdiff(rownames(OV.Beta), gene.in)
    
    strata1.in = length(intersect(gene.in, j))
    no_strata1.in = length(setdiff(gene.in, strata1.in))  
    strata1.out = length(intersect(gene.out, j))  
    no_strata1.out = length(setdiff(gene.out, strata1.out))    
    
    fish =  matrix(c(strata1.in, no_strata1.in, strata1.out, no_strata1.out), nrow = 2)    
    test = fisher.test(fish, y = NULL, alternative = "greater")  
    p.val = test$p.value     
    print(p.val)
    Metabolism[[i]][[count]] = p.val
    count = count + 1
  }
}

Metabolism.Results = as.data.frame(matrix(unlist(Metabolism), ncol = length(Metabolism.Paths)))
colnames(Metabolism.Results) = Metabolism.Paths
rownames(Metabolism.Results) = Groups.Names

## Cell-Cell
CellCell.Paths = c("Extracellular matrix organization","Type I hemidesmosome assembly","Cell junction organization","Cell-Cell communication","Apoptotic cleavage of cellular proteins","Assembly of collagen fibrils and other multimeric structures")

CellCell = list()
for(i in CellCell.Paths){
  count = 1
  for(j in Groups.List){
    gene.in = Pathways.Reactome[[i]]
    print(gene.in[1:10])
    gene.out = setdiff(rownames(OV.Beta), gene.in)
    
    strata1.in = length(intersect(gene.in, j))
    no_strata1.in = length(setdiff(gene.in, strata1.in))  
    strata1.out = length(intersect(gene.out, j))  
    no_strata1.out = length(setdiff(gene.out, strata1.out))    
    
    fish =  matrix(c(strata1.in, no_strata1.in, strata1.out, no_strata1.out), nrow = 2)    
    test = fisher.test(fish, y = NULL, alternative = "greater")  
    p.val = test$p.value     
    print(p.val)
    CellCell[[i]][[count]] = p.val
    count = count + 1
  }
}

CellCell.Results = as.data.frame(matrix(unlist(CellCell), ncol = length(CellCell.Paths)))
colnames(CellCell.Results) = CellCell.Paths
rownames(CellCell.Results) = Groups.Names

## DNA
DNA.Paths = c("Cell Cycle","Cell Cycle, Mitotic","DNA strand elongation","Gene expression (Transcription)","Synthesis of DNA","Factors involved in megakaryocyte development and platelet production")

DNA = list()
for(i in DNA.Paths){
  count = 1
  for(j in Groups.List){
    gene.in = Pathways.Reactome[[i]]
    print(gene.in[1:10])
    gene.out = setdiff(rownames(OV.Beta), gene.in)
    
    strata1.in = length(intersect(gene.in, j))
    no_strata1.in = length(setdiff(gene.in, strata1.in))  
    strata1.out = length(intersect(gene.out, j))  
    no_strata1.out = length(setdiff(gene.out, strata1.out))    
    
    fish =  matrix(c(strata1.in, no_strata1.in, strata1.out, no_strata1.out), nrow = 2)    
    test = fisher.test(fish, y = NULL, alternative = "greater")  
    p.val = test$p.value     
    print(p.val)
    DNA[[i]][[count]] = p.val
    count = count + 1
  }
}

DNA.Results = as.data.frame(matrix(unlist(DNA), ncol = length(DNA.Paths)))
colnames(DNA.Results) = DNA.Paths
rownames(DNA.Results) = Groups.Names

## Complement
Complement.Paths = c("Regulation of Complement cascade","Complement cascade","Terminal pathway of complement","Formation of Fibrin Clot (Clotting Cascade)","Intrinsic Pathway of Fibrin Clot Formation","Hemostasis","Common Pathway of Fibrin Clot Formation")

Complement = list()
for(i in Complement.Paths){
  count = 1
  for(j in Groups.List){
    gene.in = Pathways.Reactome[[i]]
    print(gene.in[1:10])
    gene.out = setdiff(rownames(OV.Beta), gene.in)
    
    strata1.in = length(intersect(gene.in, j))
    no_strata1.in = length(setdiff(gene.in, strata1.in))  
    strata1.out = length(intersect(gene.out, j))  
    no_strata1.out = length(setdiff(gene.out, strata1.out))    
    
    fish =  matrix(c(strata1.in, no_strata1.in, strata1.out, no_strata1.out), nrow = 2)    
    test = fisher.test(fish, y = NULL, alternative = "greater")  
    p.val = test$p.value     
    print(p.val)
    Complement[[i]][[count]] = p.val
    count = count + 1
  }
}

Complement.Results = as.data.frame(matrix(unlist(Complement), ncol = length(Complement.Paths)))
colnames(Complement.Results) = Complement.Paths
rownames(Complement.Results) = Groups.Names

## Cytokine
Cytokine.Paths = c("Immune System","Cytokine Signaling in Immune system","Interferon Signaling","Interferon gamma signaling","Adaptive Immune System")

Cytokine = list()
for(i in Cytokine.Paths){
  count = 1
  for(j in Groups.List){
    gene.in = Pathways.Reactome[[i]]
    gene.out = setdiff(rownames(OV.Beta), gene.in)
    
    strata1.in = length(intersect(gene.in, j))
    no_strata1.in = length(setdiff(gene.in, strata1.in))  
    strata1.out = length(intersect(gene.out, j))  
    no_strata1.out = length(setdiff(gene.out, strata1.out))    
    
    fish =  matrix(c(strata1.in, no_strata1.in, strata1.out, no_strata1.out), nrow = 2)    
    test = fisher.test(fish, y = NULL, alternative = "greater")  
    p.val = test$p.value     
    print(p.val)
    Cytokine[[i]][[count]] = p.val
    count = count + 1
  }
}

Cytokine.Results = as.data.frame(matrix(unlist(Cytokine), ncol = length(Cytokine.Paths)))
colnames(Cytokine.Results) = Cytokine.Paths
rownames(Cytokine.Results) = Groups.Names

## Erythrocyte
Erythrocyte.Paths = c("Factors involved in megakaryocyte development and platelet production","Hemostasis","Common Pathway of Fibrin Clot Formation","Interaction between L1 and Ankyrins")

Erythrocyte = list()
for(i in Erythrocyte.Paths){
  count = 1
  for(j in Groups.List){
    gene.in = Pathways.Reactome[[i]]
    gene.out = setdiff(rownames(OV.Beta), gene.in)
    
    strata1.in = length(intersect(gene.in, j))
    no_strata1.in = length(setdiff(gene.in, strata1.in))  
    strata1.out = length(intersect(gene.out, j))  
    no_strata1.out = length(setdiff(gene.out, strata1.out))    
    
    fish =  matrix(c(strata1.in, no_strata1.in, strata1.out, no_strata1.out), nrow = 2)    
    test = fisher.test(fish, y = NULL, alternative = "greater")  
    p.val = test$p.value     
    print(p.val)
    Erythrocyte[[i]][[count]] = p.val
    count = count + 1
  }
}

Erythrocyte.Results = as.data.frame(matrix(unlist(Erythrocyte), ncol = length(Erythrocyte.Paths)))
colnames(Erythrocyte.Results) = Erythrocyte.Paths
rownames(Erythrocyte.Results) = Groups.Names

## All Pathways
All.Pathway.pValues = cbind(ECM.Results, Metabolism.Results, CellCell.Results, DNA.Results, Complement.Results, Cytokine.Results, Erythrocyte.Results)
All.Pathway.pValues = as.data.frame(t(All.Pathway.pValues))
All.Pathway.pValues$Category = c("ECM interaction","ECM interaction","ECM interaction","ECM interaction","ECM interaction","ECM interaction","metabolism","metabolism","metabolism","metabolism","metabolism","cell-cell communications","cell-cell communications","cell-cell communications","cell-cell communications","cell-cell communications","cell-cell communications","DNA replication","DNA replication","DNA replication","DNA replication","DNA replication","DNA replication","complement cascade","complement cascade","complement cascade","complement cascade","complement cascade","complement cascade","complement cascade","cytokine signaling","cytokine signaling","cytokine signaling","cytokine signaling","cytokine signaling","erythrocyte and platelet","erythrocyte and platelet","erythrocyte and platelet","erythrocyte and platelet")

write.table(All.Pathway.pValues, file = "All.Pathway.pValues.Group5.txt", sep = "\t", quote = FALSE, row.names = TRUE)

############################################################
## Figure S4G
############################################################
text_list = list(text1 = "cytokine signaling",
                 text2 = "ECM interaction",
                 text3 = "DNA replication",
                 text4 = "",
                 text5 = "DNA replication\nmetabolism")
Right_Annotation_Test = rowAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(unlist(text_list)) + unit(4, "mm")))

set.seed(12345)
pdf("Pathway/FD_RNA_TCGA_DiffExp_Heatmap_AllSign_Cluster5_Labels.pdf", height = 6, width = 10, useDingbats = FALSE)
ht_list1 = Heatmap(OV.Beta.t.scale.2.t, name = "Z-score", top_annotation = FD.RNA.Ann.2, use_raster = TRUE,
                   cluster_rows = TRUE, cluster_columns = TRUE, column_split = factor(OV.Beta.meta$Verhaak),
                   border = TRUE, show_column_names = FALSE, show_row_names = FALSE, show_row_dend = FALSE,
                   row_names_gp = gpar(fontsize = 32), row_km = 5, row_title = c("","","","",""), 
                   left_annotation = c(FD.RNA.All, FD.Genes.RowAnnotation),
                   right_annotation = Right_Annotation_Test)
ht_list1 = draw(ht_list1, heatmap_legend_side = c("right"), padding = unit(c(2, 2, 2, 2), "mm"))

for(i in 1:5) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}
dev.off()

############################################################
## Figure S4H
############################################################
FD.RNA.pathway = read.table("input_pathway/FD_RNA_pathway_20220208.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)

FD.RNA.pathway.t.scale = scale(t(FD.RNA.pathway))
FD.RNA.pathway.t.scale[is.na(FD.RNA.pathway.t.scale)] <- 0
FD.RNA.pathway.t.scale.t = as.matrix(t(FD.RNA.pathway.t.scale))

FD.Figure3.meta.order.RNA = FD.Figure3.meta.order.1
FD.Figure3.meta.order.RNA = FD.Figure3.meta.order.RNA[colnames(FD.RNA.pathway),]

FD.RNA.meta.Ann = HeatmapAnnotation(sample.age = FD.Figure3.meta.order.RNA$age,
  location = FD.Figure3.meta.order.RNA$location,
  TCGA.subtype = FD.Figure3.meta.order.RNA$TCGA.subtype,
  tumor.purity = FD.Figure3.meta.order.RNA$tumor.purity,
  genome.instability = FD.Figure3.meta.order.RNA$genome.instability,
  response = FD.Figure3.meta.order.RNA$response,
  Protein.cluster = FD.Figure3.meta.order.RNA$subtype,
  RNA.cluster = FD.Figure3.meta.order.RNA$RNA_cluster,
  col = list(location = c("OV" = "yellow","OM" = "red","other" = "dodgerblue","mix" = "orange"),
    response = c("sensitive" = "gold","refractory" = "blue"),
    sample.age = FD.sample.age.col, tumor.purity = FD.tumor.purity.col,
    TCGA.subtype = c("Proliferative" = "blue","Mesenchymal" = "violet","Immunoreactive" = "red","Differentiated" = "green"),
    Protein.cluster = c("1" = "#F8766D","2" = "#BB9D00","3" = "#00B81F","4" = "#00A5FF","5" = "#E76BF3"),
    RNA.cluster = c("1" = "#F8766D","2" = "#BB9D00","3" = "#00B81F","4" = "#00A5FF","5" = "#E76BF3")), na_col = "lightgray")

theme_set(theme_grey(base_size = 16))
pdf("Figure_3/FD_RNA_Pathway_Test.pdf", height = 14, width = 16, useDingbats = FALSE)
ht_list = Heatmap(FD.RNA.pathway.t.scale.t, top_annotation = FD.RNA.meta.Ann,
                  cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(FD.Figure3.meta.order.RNA$RNA_cluster),
                  row_split = factor(FD.Figure3.pathwayGroups$Group, levels = c("6","5","4","3","2","1")),
                  border = TRUE, 
                  row_title = c("","","","","",""), 
                  column_title = c("","","","",""),
                  show_column_names = FALSE, row_labels = FD.label.pathways.2, left_annotation = FD.Figure3.meta.Row)
draw(ht_list, heatmap_legend_side = c("bottom"), annotation_legend_side = "bottom",  padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()
