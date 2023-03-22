############################################################
## title: "Figure_2_S2"
## author: "Oscar Murillo"
## date: "August 18, 2020"
## panels: "2G, 2H, 2I, 2J, S2D, S2E, S2F, S2G, S2H, S2J"
############################################################
FD.protein.df = FD.protein
FD.protein.df = as.data.frame(t(FD.protein.df))
FD.protein.df$Sample = rownames(FD.protein.df)
FD.protein.df.TP53 = FD.protein.df[,c("Sample","TP53","CDC20","KIF2C","PLK1")]
rownames(FD.protein.df.TP53) = NULL
colnames(FD.protein.df.TP53) = c("Sample","TP53_protein","CDC20_protein","KIF2C_protein","PLK1_protein")

FD.RNA = read.table("input_RNA/RNAcounts_normalized_FFPE_2021521.txt", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = 1)
FD.RNA.df = FD.RNA
FD.RNA.df = as.data.frame(t(FD.RNA.df))
FD.RNA.df$Sample = rownames(FD.RNA.df)
FD.RNA.df.TP53 = FD.RNA.df[,c("Sample","TP53","MDM2","CDKN1A","CDC20","CENPA","KIF2C")]
rownames(FD.RNA.df.TP53) = NULL
colnames(FD.RNA.df.TP53) = c("Sample","TP53_RNA","MDM2_RNA","CDKN1A_RNA","CDC20_RNA","CENPA_RNA","KIF2C_RNA")

FD.LOH.data = read.table("LOH.txt", header = TRUE, sep = "\t", quote = NULL, row.names = NULL)
FD.LOH.data$chr17LOH = factor(FD.LOH.data$chr17LOH)
FD.LOH.data = FD.LOH.data[,c("Sample","chr17LOH")]

FD.mutations.data = read.table("input_TP53/FD_DNA_mutations.txt", header = TRUE, sep = "\t", quote = NULL, row.names = NULL)
FD.mutations.data = FD.mutations.data[,c("Sample","Alterations","IMPACT")]
FD.mutations.data$Alterations = gsub("frameshift_variant","Truncating", FD.mutations.data$Alterations)
FD.mutations.data$Alterations = gsub("inframe_deletion","Missense", FD.mutations.data$Alterations)
FD.mutations.data$Alterations = gsub("splice_acceptor_variant","Truncating", FD.mutations.data$Alterations)
FD.mutations.data$Alterations = gsub("splice_donor_variant","Truncating", FD.mutations.data$Alterations)
FD.mutations.data$Alterations = gsub("missense_variant","Missense", FD.mutations.data$Alterations)
FD.mutations.data$Alterations = gsub("stop_gained","Truncating", FD.mutations.data$Alterations)

FD.Combined.1 = merge(FD.meta, FD.protein.df.TP53, by = "Sample", all.x = TRUE)
FD.Combined.2 = merge(FD.Combined.1, FD.RNA.df.TP53, by = "Sample", all.x = TRUE)
FD.Combined.3 = FD.Combined.2
FD.Combined.4 = merge(FD.Combined.3, FD.LOH.data, by = "Sample", all.x = TRUE)
FD.Combined.5 = merge(FD.Combined.4, FD.mutations.data, by = "Sample", all.x = TRUE)
FD.Combined.5$Alterations[is.na(FD.Combined.5$Alterations)] <- "WT"
FD.Combined.5$Alterations = factor(FD.Combined.5$Alterations, levels = c("WT","Missense","Truncating"))
FD.Combined.5$IMPACT[is.na(FD.Combined.5$IMPACT)] <- "WT"
FD.Combined.5$IMPACT = factor(FD.Combined.5$IMPACT, levels = c("WT","MODERATE","HIGH"))

FD.Combined.5 = FD.Combined.5 %>% group_by(Sample) %>% filter(duplicated(Sample) | n() == 1)
FD.Combined.5.melt = melt(FD.Combined.5)

## Figure S2D
FD.TP53.protein.mutations = ggplot(data = FD.Combined.5 %>% filter(!is.na(chr17LOH)), aes(x = Alterations, y = TP53_protein)) + 
  geom_boxplot(outlier.shape = NA) + xlab("TP53 alterations") + ylab("TP53 protein") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", ref.group = "WT") +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

## Figure S2G-1
FD.MDM2.RNA.LOH = ggplot(data = FD.Combined.5 %>% filter(!is.na(chr17LOH)), aes(x = chr17LOH, y = MDM2_RNA)) + 
  geom_boxplot(outlier.shape = NA) + xlab("chr17LOH") + ylab("MDM2 RNA") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", label.x = 1.25) +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

## Figure S2G-2
FD.CDC20_RNA.LOH = ggplot(data = FD.Combined.5 %>% filter(!is.na(chr17LOH)), aes(x = chr17LOH, y = CDC20_RNA)) + 
  geom_boxplot(outlier.shape = NA) + xlab("chr17LOH") + ylab("CDC20 RNA") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", label.x = 1.25) +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

## Figure 2G
FD.TP53.Protein.LOH = ggplot(data = FD.Combined.5 %>% filter(!is.na(chr17LOH)), aes(x = chr17LOH, y = TP53_protein)) + 
  geom_boxplot(outlier.shape = NA) + xlab("chr17LOH") + ylab("TP53 protein") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", label.x = 1.25) +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

###############################################
## Activity Scores
###############################################
FD.Combined.Test = FD.Combined.5
FD.Combined.Test = FD.Combined.Test[,c("Sample","Response","Type","Location","Age","Neo","Stage","Source","Sample_Age","Cluster","TP53_protein","TP53_RNA","TP53_ActivityScore","chr17LOH","Alterations","IMPACT" )]
FD.RNA.Combined.Test = FD.Combined.Test[FD.Combined.Test$Sample %in% colnames(FD.RNA),]

TP53.WT.Sign = c("ABCA12","APOBEC3H","ASCC3","BAX","BBC3","BTG2","C9orf169","CCNG1","CDKN1A","DCP1B","FAS","GADD45A","GPR87","KIAA0247","LIF","MDM2","NUPR1","PGF","PLK2","PPM1D","PSTPIP2","PTP4A1","RAP2B","RPS27L","SESN1","TNFRSF10B","TRAF4","TRIAP1","TSKU","ZMAT3","ZNF79")
TP53.Mutate.Sign = c("AUNIP","BUB1","CDC20","CDCA5","CDCA8","CENPA","CENPI","DDIAS","DEPDC1","ERCC6L","FAM72B","KIF2C","MELK","NDC80","OIP5","POLQ","PLK1","SGOL1","TPX2","TTK")

###############################################
## ssGSEA Scores: RNA
###############################################
TP53.Genes = list(TP53.WT.Sign, TP53.Mutate.Sign)
names(TP53.Genes) = c("ssGSEA.TP53.WT","ssGSEA.TP53.Mut")

FD.RNA.TP53.GSEA = gsva(expr = data.matrix(FD.RNA), TP53.Genes, method = "ssgsea")
FD.RNA.TP53.GSEA = as.data.frame(FD.RNA.TP53.GSEA)
FD.RNA.TP53.GSEA = as.data.frame(t(as.data.frame(FD.RNA.TP53.GSEA)))
FD.RNA.TP53.GSEA$Sample = rownames(FD.RNA.TP53.GSEA)
colnames(FD.RNA.TP53.GSEA) = c("ssGSEA.TP53.WT","ssGSEA.TP53.Mut","Sample")

## Scores are used in the analysis described in Figure 2J
FD.RNA.Combined.Test = FD.Combined.Test[FD.Combined.Test$Sample %in% colnames(FD.RNA),]
FD.RNA.Combined.Test.GSEA = merge(FD.RNA.Combined.Test, FD.RNA.TP53.GSEA, by = "Sample")

## Figure 2H
FD.RNA.GSEA.WT.LOH = ggplot(data = FD.RNA.Combined.Test.GSEA %>% filter(!is.na(chr17LOH)), aes(x = chr17LOH, y = ssGSEA.TP53.WT)) + 
  geom_boxplot(outlier.shape = NA) + xlab("chr17LOH") + ylab("WT TP53 score") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

## Figure 2I
FD.RNA.GSEA.WT.Response = ggplot(data = FD.RNA.Combined.Test.GSEA, aes(x = Response, y = ssGSEA.TP53.WT, color = Response)) + 
  geom_boxplot(outlier.shape = NA) + xlab("Response") + ylab("WT TP53 score") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", label = "p.format") + scale_color_manual(values = c("#0073C2","#EFC000")) +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

## Figure S2F
FD.RNA.GSEA.Mutated.LOH = ggplot(data = FD.RNA.Combined.Test.GSEA %>% filter(!is.na(chr17LOH)), aes(x = chr17LOH, y = ssGSEA.TP53.Mut)) + 
  geom_boxplot(outlier.shape = NA) + xlab("chr17-LOH") + ylab("TP53-mutant activity score (RNA)") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

## Figure S2H
FD.RNA.GSEA.Mutated.Response = ggplot(data = FD.RNA.Combined.Test.GSEA, aes(x = Response, y = ssGSEA.TP53.Mut, color = Response)) + 
  geom_boxplot(outlier.shape = NA) + xlab("response") + ylab("TP53-mutant activity score (RNA)") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  stat_compare_means(method = "wilcox.test", label = "p.format") + scale_color_manual(values = c("#0073C2","#EFC000")) +
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

FD.RNA.Combined.Test.GSEA$Alterations_Cat = FD.RNA.Combined.Test.GSEA$Alterations
FD.RNA.Combined.Test.GSEA$Alterations_Cat = gsub("Missense", "Mutated", FD.RNA.Combined.Test.GSEA$Alterations_Cat)
FD.RNA.Combined.Test.GSEA$Alterations_Cat = gsub("Truncating", "Mutated", FD.RNA.Combined.Test.GSEA$Alterations_Cat)

## Figure S2E
FD.RNA.GSEA.WT.Alterations = ggplot(data = FD.RNA.Combined.Test.GSEA %>% filter(!is.na(chr17LOH)), aes(x = Alterations_Cat, y = ssGSEA.TP53.WT)) + 
  geom_boxplot(outlier.shape = NA) + xlab("TP status") + ylab("TP53-WT activity score (RNA)") + theme_bw() + theme(legend.position = "none") + geom_jitter(shape = 16, position = position_jitter(0.2)) + 
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  scale_y_continuous(labels = scaleFUN)

###############################################
## Figure S2J
###############################################
FD.RNA.Combined.Test.GSEA$chr17LOH = factor(FD.RNA.Combined.Test.GSEA$chr17LOH, levels = c(1,0))

FD.BRCA = read.table("input_TP53/BRCA_Mutations.txt", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = NULL)

FD.RNA.Combined.Test.GSEA.BRCA = merge(FD.RNA.Combined.Test.GSEA, FD.BRCA, by = "Sample", all.x = TRUE)
FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH = factor(FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH)
FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH = ifelse(FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH == "0", yes = "2", no = "1")
FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH[is.na(FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH)] <- 0
FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH = factor(FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH)

FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH = gsub("0","Unknown", FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH)
FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH = gsub("1","Yes", FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH)
FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH = gsub("2","No", FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH)
FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH = factor(FD.RNA.Combined.Test.GSEA.BRCA$chr17LOH, levels = c("Yes","No","Unknown"))

FD.RNA.GSEA.XY = ggplot(FD.RNA.Combined.Test.GSEA.BRCA, aes(x = ssGSEA.TP53.WT, y = ssGSEA.TP53.Mut, color = Response, group = Response, shape = chr17LOH, label = Gene)) +
  labs(x = "WT TP53 score", y = "Mutated TP53 score") + theme_bw() + scale_color_manual(values = c("#0073C2","#EFC000")) + 
  theme(text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  geom_point(size = 3) + stat_cor(label.x = 0.1, label.y = c(0.4, 0.45), size = 4.5) +
  geom_smooth(method = lm, fullrange = FALSE, formula = y~x, se = FALSE) + theme(legend.position = "bottom", plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_text(hjust = 0.5, nudge_x = 0, nudge_y = 0.03, size = 3)
