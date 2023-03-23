############################################################
## title: "Figure 5A_5C_5D_S5A_S5D"
## author: "Oscar Murillo"
## date: "March 22, 2023"
## panels: "5A, 5C, 5D, S5A, S5D"
############################################################
############################################################
## Read in pathway datasets
############################################################
FD.pathway = read.table("input_pathway/FD_pathway_20210927.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)
FZ.pathway = read.table("input_pathway/FZ_pathway_20210927.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)
Retro.pathway = read.table("input_pathway/CPTAC2_retro_pathway_20211004.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)
PDX.pathway = read.table("input_pathway/PDX_pathway_20211216.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)

FD.meta = read.table("input_metadata/FD_GLBLprot_metadata.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)
FZ.meta = read.table("input_metadata/FZ_GLBLprot_metadata.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)
Retro.meta = read.table("input_metadata/CPTAC2_retro_metadata.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)

############################################################
## Perform PAM training using the FD clusters
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
pamr.confusion(FD.pathway.results, threshold = FD.pathway.threshold)

############################################################
## Determine FD cluster centroids for heatmap
############################################################
FD.pathway.top = pamr.listgenes(FD.pathway.training, FD.pathway.input, threshold = FD.pathway.threshold)
FD.pathway.top = as.data.frame(FD.pathway.top)
rownames(FD.pathway.top) = FD.pathway.top$id
FD.pathway.top$id = NULL
colnames(FD.pathway.top) = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")
FD.pathway.top$Pathway = rownames(FD.pathway.top)
FD.pathway.top = FD.pathway.top[rev(pathway.order),]
pathway.order = rownames(FD.pathway)

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

############################################################
## Predict PAM classification for FZ
############################################################
FD.pathway.temp1 = FD.pathway
FD.pathway.temp1$Mean = rowMeans(FD.pathway.temp1)
FD.pathway.temp1$SD = apply(FD.pathway.temp1[,1:length(colnames(FD.pathway))], 1, sd)

FZ.pathway.temp1 = FZ.pathway
FZ.pathway.temp1$Mean = rowMeans(FZ.pathway.temp1)
FZ.pathway.temp1$SD = apply(FZ.pathway.temp1[,1:length(colnames(FZ.pathway))], 1, sd)

FZ.pathway.norm = NULL
for(i in 1:150){
  FZ.pathway.norm[[i]] = ((FZ.pathway[i,] - FZ.pathway.temp1$Mean[i])*(FD.pathway.temp1$SD[i]/FZ.pathway.temp1$SD[i])) + FD.pathway.temp1$Mean[i]
}

FZ.pathway.norm.df = FZ.pathway.norm
FZ.pathway.norm.df = data.frame(matrix(unlist(FZ.pathway.norm.df), nrow = length(FZ.pathway.norm.df), byrow = TRUE))
colnames(FZ.pathway.norm.df) = colnames(FZ.pathway)
rownames(FZ.pathway.norm.df) = rownames(FZ.pathway)

FD.pathway.overlap = FD.pathway[,FD.sampleOverlap]
FD.pathway.overlap$Mean = rowMeans(FD.pathway.overlap)
FD.pathway.overlap$SD = apply(FD.pathway.overlap[,1:length(FD.sampleOverlap)], 1, sd)
FZ.pathway.overlap = FZ.pathway[,FZ.sampleOverlap]
FZ.pathway.overlap$Mean = rowMeans(FZ.pathway.overlap)
FZ.pathway.overlap$SD = apply(FZ.pathway.overlap[,1:length(FZ.sampleOverlap)], 1, sd)

FZ.pathway.overlap.norm = NULL
for(i in 1:150){
  FZ.pathway.overlap.norm[[i]] = ((FZ.pathway[i,] - FZ.pathway.overlap$Mean[i])*(FD.pathway.overlap$SD[i]/FZ.pathway.overlap$SD[i])) + FD.pathway.overlap$Mean[i]
}

FZ.pathway.overlap.norm.df = FZ.pathway.overlap.norm
FZ.pathway.overlap.norm.df = data.frame(matrix(unlist(FZ.pathway.overlap.norm.df), nrow = length(FZ.pathway.overlap.norm.df), byrow = TRUE))
colnames(FZ.pathway.overlap.norm.df) = colnames(FZ.pathway)
rownames(FZ.pathway.overlap.norm.df) = rownames(FZ.pathway)

FZ.pathway.matrix = as.matrix(FZ.pathway.overlap.norm.df)
FZ.pathway.y = NULL
FZ.pathway.genenames = rownames(FZ.pathway.overlap.norm.df)
FZ.pathway.geneid = rownames(FZ.pathway.overlap.norm.df)
FZ.pathway.samples = colnames(FZ.pathway.overlap.norm.df)

FZ.pathway.input = list(x = FZ.pathway.matrix, y = FZ.pathway.y, genenames = FZ.pathway.genenames, geneid = FZ.pathway.geneid, sampleid = FZ.pathway.samples)
FZ.pathway.predict = pamr.predict(FD.pathway.training, FZ.pathway.input$x, threshold = FD.pathway.threshold, type = "posterior")
FZ.pathway.predict = as.data.frame(FZ.pathway.predict)
FZ.pathway.predict$PAM = colnames(FZ.pathway.predict)[apply(FZ.pathway.predict, 1, which.max)]
FZ.pathway.predict$Sample = FZ.pathway.samples
FZ.pathway.predict = FZ.pathway.predict %>% dplyr::select(Sample, everything())

write.table(FZ.pathway.predict, file = "PAM/FZ_Pathway_Probability_V1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Predict PAM classification for CPTAC
############################################################
FD.pathway.temp1 = FD.pathway
FD.pathway.temp1$Mean = rowMeans(FD.pathway.temp1)
FD.pathway.temp1$SD = apply(FD.pathway.temp1[,1:length(colnames(FD.pathway))], 1, sd)

Retro.pathway.temp1 = Retro.pathway
Retro.pathway.temp1$Mean = rowMeans(Retro.pathway.temp1)
Retro.pathway.temp1$SD = apply(Retro.pathway.temp1[,1:length(colnames(Retro.pathway))], 1, sd)

Retro.pathway.norm = NULL
for(i in 1:150){
  Retro.pathway.norm[[i]] = ((Retro.pathway[i,] - Retro.pathway.temp1$Mean[i])*(FD.pathway.temp1$SD[i]/Retro.pathway.temp1$SD[i])) + FD.pathway.temp1$Mean[i]
}

Retro.pathway.norm.df = Retro.pathway.norm
Retro.pathway.norm.df = data.frame(matrix(unlist(Retro.pathway.norm.df), nrow = length(Retro.pathway.norm.df), byrow = TRUE))
colnames(Retro.pathway.norm.df) = colnames(Retro.pathway)
rownames(Retro.pathway.norm.df) = rownames(Retro.pathway)

Retro.pathway.matrix = as.matrix(Retro.pathway.norm.df)
Retro.pathway.y = NULL
Retro.pathway.genenames = rownames(Retro.pathway.norm.df)
Retro.pathway.geneid = rownames(Retro.pathway.norm.df)
Retro.pathway.samples = colnames(Retro.pathway.norm.df)

Retro.pathway.input = list(x = Retro.pathway.matrix, y = Retro.pathway.y, genenames = Retro.pathway.genenames, geneid = Retro.pathway.geneid, sampleid = Retro.pathway.samples)
Retro.pathway.predict = pamr.predict(FD.pathway.training, Retro.pathway.input$x, threshold = FD.pathway.threshold, type = "posterior")
Retro.pathway.predict = as.data.frame(Retro.pathway.predict)
Retro.pathway.predict$PAM = colnames(Retro.pathway.predict)[apply(Retro.pathway.predict, 1, which.max)]
Retro.pathway.predict$Sample = Retro.pathway.samples
Retro.pathway.predict = Retro.pathway.predict %>% dplyr::select(Sample, everything())

write.table(Retro.pathway.predict, file = "PAM/Retro_Pathway_Probability.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Predict PAM classification for PDX 
############################################################
FD.pathway.temp4 = FD.pathway
FD.pathway.temp4$Mean = rowMeans(FD.pathway.temp4)
FD.pathway.temp4$SD = apply(FD.pathway.temp4[,1:length(colnames(FD.pathway))], 1, sd)

PDX.pathway.temp4 = PDX.pathway
PDX.pathway.temp4$Mean = rowMeans(PDX.pathway.temp4)
PDX.pathway.temp4$SD = apply(PDX.pathway.temp4[,1:length(colnames(PDX.pathway))], 1, sd)

PDX.pathway.norm = NULL
for(i in 1:150){
  PDX.pathway.norm[[i]] = ((PDX.pathway[i,] - PDX.pathway.temp4$Mean[i])*(FD.pathway.temp4$SD[i]/PDX.pathway.temp4$SD[i])) + FD.pathway.temp4$Mean[i]
}

PDX.pathway.norm.df = PDX.pathway.norm
PDX.pathway.norm.df = data.frame(matrix(unlist(PDX.pathway.norm.df), nrow = length(PDX.pathway.norm.df), byrow = TRUE))
colnames(PDX.pathway.norm.df) = colnames(PDX.pathway)
rownames(PDX.pathway.norm.df) = rownames(PDX.pathway)
PDX.pathway.norm.df = PDX.pathway.norm.df[, order(colSums(PDX.pathway.norm.df), decreasing = TRUE)]

PDX.pathway.matrix = as.matrix(PDX.pathway.norm.df)
PDX.pathway.y = NULL
PDX.pathway.genenames = rownames(PDX.pathway.norm.df)
PDX.pathway.geneid = rownames(PDX.pathway.norm.df)
PDX.pathway.samples = colnames(PDX.pathway.norm.df)

PDX.pathway.input = list(x = PDX.pathway.matrix, y = PDX.pathway.y, genenames = PDX.pathway.genenames, geneid = PDX.pathway.geneid, sampleid = PDX.pathway.samples)
PDX.pathway.predict = pamr.predict(FD.pathway.training, PDX.pathway.input$x, threshold = FD.pathway.threshold, type = "posterior")
PDX.pathway.predict = as.data.frame(PDX.pathway.predict)
PDX.pathway.predict$PAM = colnames(PDX.pathway.predict)[apply(PDX.pathway.predict, 1, which.max)]
PDX.pathway.predict$Sample = PDX.pathway.samples
PDX.pathway.predict = PDX.pathway.predict %>% dplyr::select(Sample, everything())

write.table(PDX.pathway.predict, file = "PAM/PDX_Pathway_Probability.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Structure FZ (non-overlap) for combined heatmap
############################################################
FZ.pathway.order = FZ.pathway.overlap.norm.df[, FZ.pathway.predict.order$Sample]
FZ.meta.order = FZ.meta
rownames(FZ.meta.order) = FZ.meta.order$Sample
FZ.meta.order = FZ.meta.order[FZ.pathway.predict.order$Sample,]

FZ.pathway.order.t.scale = scale(t(FZ.pathway.order))
FZ.pathway.order.t.scale[is.na(FZ.pathway.order.t.scale)] <- 0
FZ.pathway.order.t.scale.t = as.matrix(t(FZ.pathway.order.t.scale))

FZ.Figure3.meta = read.table("Figure_3/FZ_metadata.txt", sep = "\t", header = TRUE)
rownames(FZ.Figure3.meta) = FZ.Figure3.meta$Sample
FZ.Figure3.meta = FZ.Figure3.meta[FZ.pathway.predict.order$Sample,]
FZ.Figure3.meta$subtype = factor(FZ.Figure3.meta$subtype, levels = c(1,2,3,4,5))

FZ.Hierc = read.table("Figure_3/km7Frozen.k=5.consensusClass.txt", sep = "\t", header = TRUE)
FZ.Hierc$hierarchical.clustering = as.factor(FZ.Hierc$hierarchical.clustering)
FZ.Figure3.meta.Hierc = merge(FZ.Figure3.meta, FZ.Hierc, by = "Sample")
rownames(FZ.Figure3.meta.Hierc) = FZ.Figure3.meta.Hierc$Sample
FZ.Figure3.meta.Hierc = FZ.Figure3.meta.Hierc[colnames(FZ.pathway.order.t.scale.t),]
FZ.Figure3.meta.Hierc$overlap29 = ifelse(FZ.Figure3.meta.Hierc$Sample %in% FZ.sampleOverlap, "Overlap", "No_Overlap")
FZ.Figure3.meta.Hierc.NotOverlap = FZ.Figure3.meta.Hierc[FZ.Figure3.meta.Hierc$overlap29 == "No_Overlap",]
FZ.Figure3.meta.Hierc.NotOverlap = FZ.Figure3.meta.Hierc.NotOverlap[order(FZ.Figure3.meta.Hierc.NotOverlap[,7], FZ.Figure3.meta.Hierc.NotOverlap[,11], FZ.Figure3.meta.Hierc.NotOverlap[,2]),]

FZ.pathway.order.t.scale.t.NotOverlap = FZ.pathway.order.t.scale.t[,FZ.Figure3.meta.Hierc.NotOverlap$Sample]

############################################################
## Structure CPTAC for combined heatmap
############################################################
Retro.pathway.order = Retro.pathway.norm.df[, Retro.pathway.predict.order$Sample]
Retro.meta.order = Retro.meta
rownames(Retro.meta.order) = Retro.meta.order$Sample
Retro.meta.order = Retro.meta.order[Retro.pathway.predict.order$Sample,]

Retro.pathway.order.t.scale = scale(t(Retro.pathway.order))
Retro.pathway.order.t.scale[is.na(Retro.pathway.order.t.scale)] <- 0
Retro.pathway.order.t.scale.t = as.matrix(t(Retro.pathway.order.t.scale))

Retro.Figure3.meta.order = Retro.pathway.predict.order
Retro.Figure3.meta.order$PAM = gsub("Cluster_","",Retro.Figure3.meta.order$PAM)
Retro.Figure3.meta.order$PAM = factor(Retro.Figure3.meta.order$PAM, levels = c(1,2,3,4,5))

Retro.Hierc = read.table("Figure_3/Retro_hiera.txt", sep = "\t", header = TRUE)
Retro.Figure3.meta.order.2 = merge(Retro.Figure3.meta.order, Retro.Hierc, by = "Sample")
Retro.Figure3.meta.order.2$hierarchical.clustering = factor(Retro.Figure3.meta.order.2$hierarchical.clustering, levels = c(1,2,3,4,5))
rownames(Retro.Figure3.meta.order.2) = Retro.Figure3.meta.order.2$Sample
Retro.Figure3.meta.order.2 = Retro.Figure3.meta.order.2[colnames(Retro.pathway.order.t.scale.t),]

Retro.Figure3.meta.order.2 = Retro.Figure3.meta.order.2[order(Retro.Figure3.meta.order.2[,7], Retro.Figure3.meta.order.2[,9]),]
Retro.pathway.order.t.scale.t = Retro.pathway.order.t.scale.t[,Retro.Figure3.meta.order.2$Sample]

############################################################
## Structure PDX for combined heatmap
############################################################
PDX.pathway.order = PDX.pathway.norm.df[, PDX.pathway.predict.order$Sample]
PDX.pathway.order.t.scale = scale(t(PDX.pathway.order))
PDX.pathway.order.t.scale[is.na(PDX.pathway.order.t.scale)] <- 0
PDX.pathway.order.t.scale.t = as.matrix(t(PDX.pathway.order.t.scale))

PDX.Figure3.meta.order = PDX.pathway.predict.order
PDX.Figure3.meta.order$PAM = gsub("Cluster_","",PDX.Figure3.meta.order$PAM)
PDX.Figure3.meta.order$PAM = factor(PDX.Figure3.meta.order$PAM, levels = c(1,2,3,4,5))

PDX.Hierc = read.table("Figure_3/km7PDX.k=5.consensusClass.txt", sep = "\t", header = TRUE)
PDX.Figure3.meta.order.2 = merge(PDX.Figure3.meta.order, PDX.Hierc, by = "Sample")
PDX.Figure3.meta.order.2$hierarchical.clustering = factor(PDX.Figure3.meta.order.2$hierarchical.clustering, levels = c(1,2,3,4,5))
rownames(PDX.Figure3.meta.order.2) = PDX.Figure3.meta.order.2$Sample
PDX.Figure3.meta.order.2 = PDX.Figure3.meta.order.2[colnames(PDX.pathway.order.t.scale.t),]

############################################################
## Figure 5A
############################################################
FZ.Figure3.meta.Ann.NonOver = HeatmapAnnotation(clustering = FZ.Figure3.meta.Hierc.NotOverlap$hierarchical.clustering, subtype = FZ.Figure3.meta.Hierc.NotOverlap$subtype,
                                                col = list(clustering = c("1" = "#F8766D","2" = "#00B81F","3" = "#BB9D00","4" = "#00A5FF","5" = "#E76BF3"),
                                                           subtype = c("1" = "#F8766D","2" = "#BB9D00","3" = "#00B81F","4" = "#00A5FF","5" = "#E76BF3")),
                                                na_col = "lightgray", annotation_name_side = "left")

Retro.Figure3.meta.Ann = HeatmapAnnotation(clustering = Retro.Figure3.meta.order.2$hierarchical.clustering, subtype = Retro.Figure3.meta.order.2$PAM,
                                           col = list(clustering = c("1" = "#00B81F","2" = "#00A5FF","3" = "#F8766D","4" = "#E76BF3","5" = "#BB9D00"),
                                                      subtype = c("1" = "#F8766D","2" = "#BB9D00","3" = "#00B81F","4" = "#00A5FF","5" = "#E76BF3")),
                                           na_col = "lightgray", show_legend = FALSE, show_annotation_name = c(clustering = FALSE, subtype = FALSE))

PDX.Figure3.meta.Ann = HeatmapAnnotation(clustering = PDX.Figure3.meta.order.2$hierarchical.clustering, subtype = PDX.Figure3.meta.order.2$PAM,
                                         col = list(clustering = c("1" = "#E76BF3","2" = "#BB9D00","3" = "#F8766D","4" = "#00B81F","5" = "#00A5FF"),
                                                    subtype = c("1" = "#F8766D","2" = "#BB9D00","3" = "#00B81F","4" = "#00A5FF","5" = "#E76BF3")),
                                         na_col = "lightgray",show_legend = FALSE)

FZ_ht_list = Heatmap(FZ.pathway.order.t.scale.t.NotOverlap, top_annotation = FZ.Figure3.meta.Ann.NonOver,
                     cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(FZ.Figure3.meta.Hierc.NotOverlap$hierarchical.clustering, levels = c("1","3","2","4","5")),
                     row_split = factor(FD.Figure3.pathwayGroups$Group, levels = c("6","5","4","3","2","1")), border = TRUE, 
                     row_title = c("","","","","",""), column_title = c("","","","",""),
                     show_column_names = FALSE, show_row_names = FALSE, width = unit(5, "cm"))

Retro_ht_list = Heatmap(Retro.pathway.order.t.scale.t, top_annotation = Retro.Figure3.meta.Ann,
                        cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(Retro.Figure3.meta.order.2$hierarchical.clustering, levels = c("3","5","1","2","4")),
                        row_split = factor(FD.Figure3.pathwayGroups$Group, levels = c("6","5","4","3","2","1")), border = TRUE, 
                        row_title = c("","","","","",""), column_title = c("","","","",""),
                        show_column_names = FALSE, show_row_names = FALSE, left_annotation = FD.Figure3.meta.Row, width = unit(17.4, "cm"))

PDX_ht_list = Heatmap(PDX.pathway.order.t.scale.t, top_annotation = PDX.Figure3.meta.Ann,
                      cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(PDX.Figure3.meta.order.2$hierarchical.clustering, levels = c("3","2","4","5","1")),
                      row_split = factor(FD.Figure3.pathwayGroups$Group, levels = c("6","5","4","3","2","1")), border = TRUE, 
                      row_title = c("","","","","",""), column_title = c("","","","",""),
                      show_column_names = FALSE, row_labels = FD.label.pathways.2, width = unit(6, "cm"))

theme_set(theme_grey(base_size = 20))
pdf("Combine_FZ_Retro_PDX_20220919.pdf", height = 12, width = 20, useDingbats = FALSE)
ht_list_combo = Retro_ht_list + FZ_ht_list + PDX_ht_list
draw(ht_list_combo, ht_gap = unit(0.25, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom",  padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

############################################################
## Figure 5C
############################################################
FD.cluster.pca = as.data.frame(t(FD.pathway.order))
FD.cluster.pch = replicate(158, 21)
FD.order.clusters = FD.meta.order$Cluster
FD.order.clusters.1 = rep("white", length(FD.order.clusters))
FD.order.clusters.1[FD.order.clusters=="Cluster_1"] = "#F8766D"
FD.order.clusters.1[FD.order.clusters=="Cluster_2"] = "#A3A500"
FD.order.clusters.1[FD.order.clusters=="Cluster_3"] = "#00BD7A"
FD.order.clusters.1[FD.order.clusters=="Cluster_4"] = "#00B0F6"
FD.order.clusters.1[FD.order.clusters=="Cluster_5"] = "#E76BF3"

Retro.newRows = Retro.Figure3.meta.order.2[,c(1,7,9)]
Retro.vec1 = Retro.newRows$hierarchical.clustering
Retro.vec1 = gsub("1","Cluster_1",Retro.vec1) 
Retro.vec1 = gsub("2","Cluster_2",Retro.vec1) 
Retro.vec1 = gsub("3","Cluster_3",Retro.vec1) 
Retro.vec1 = gsub("4","Cluster_4",Retro.vec1) 
Retro.vec1 = gsub("5","Cluster_5",Retro.vec1) 

Retro.vec1 = gsub("Cluster_3","#F8766D",Retro.vec1) 
Retro.vec1 = gsub("Cluster_5","#A3A500",Retro.vec1) 
Retro.vec1 = gsub("Cluster_1","#00BD7A",Retro.vec1) 
Retro.vec1 = gsub("Cluster_2","#00B0F6",Retro.vec1) 
Retro.vec1 = gsub("Cluster_4","#E76BF3",Retro.vec1) 

Retro.cluster.pca = as.data.frame(t(Retro.pathway.order))
Retro.cluster.pca.new = Retro.cluster.pca[Retro.newRows$Sample,]
Retro.cluster.pch = replicate(174, 23)

FD.Retro.cluster.tSNE = rbind(FD.cluster.pca, Retro.cluster.pca.new)
FD.Retro.cluster.pch = c(FD.cluster.pch, Retro.cluster.pch)
FD.Retro.cluster.col = c(FD.order.clusters.1, Retro.vec1)

FD.Retro.cluster.tSNE.tSNE = Rtsne(FD.Retro.cluster.tSNE, perplexity = 25, pca_scale = TRUE, check_duplicates = FALSE)

theme_set(theme_grey(base_size = 20))
pdf("tSNE_FD_Retro_Pathway.pdf", height = 8, width = 8)
plot(FD.Retro.cluster.tSNE.tSNE$Y, col = FD.Retro.cluster.col, bg = FD.Retro.cluster.col,
     pch = FD.Retro.cluster.pch, cex = 1.2)
dev.off()

############################################################
## Figure 5D
############################################################
## FD Boxplots
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

ribosome.path = c("RPL11","RPL23","RPL30","RPL36","RPS3")
EMT.path = c("BGN","COL6A2","FN1","LGALS1","MMP2")
Two.paths = c(ribosome.path, EMT.path)

theme_set(theme_grey(base_size = 20))
pdf("FD_Ribosome.pdf", width = 8, height = 4)
ggboxplot(FD.protein.all.order.melt[FD.protein.all.order.melt$variable %in% ribosome.path,], x = "variable", y = "value", fill = "subtype",
          ylab = "relative protein abundance", xlab = "", palette = c("#F8766D","#A3A500","#00BD7A","#00B0F6","#E76BF3"), short.panel.labs = TRUE) +
  theme(text = element_text(size = 16))
dev.off()

theme_set(theme_grey(base_size = 20))
pdf("FD_EMT.pdf", width = 8, height = 4)
ggboxplot(FD.protein.all.order.melt[FD.protein.all.order.melt$variable %in% EMT.path,], x = "variable", y = "value", fill = "subtype",
          ylab = "relative protein abundance", xlab = "", palette = c("#F8766D","#A3A500","#00BD7A","#00B0F6","#E76BF3"), short.panel.labs = TRUE) +
  theme(text = element_text(size = 16))
dev.off()

## CPTAC Boxplots
Retro.protein = read.table("input_protein/CPTAC2_retro_alignedToFD_09072021_Ave.tsv", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = 1)
Retro.protein.all.order = Retro.protein[Two.paths, Retro.pathway.predict.order$Sample]
Retro.protein.all.order.t.scale = scale(t(Retro.protein.all.order))
Retro.protein.all.order.t.scale[is.na(Retro.protein.all.order.t.scale)] <- 0
Retro.protein.all.order.t.scale.t = as.matrix(t(Retro.protein.all.order.t.scale))
Retro.protein.all.order = as.data.frame(t(Retro.protein.all.order.t.scale.t))
Retro.protein.all.order$subtype = Retro.pathway.predict.order$PAM
Retro.protein.all.order$subtype = gsub("Cluster_","",Retro.protein.all.order$subtype)
Retro.protein.all.order$subtype = factor(Retro.protein.all.order$subtype, levels = c(1,2,3,4,5))
Retro.protein.all.order$Sample = rownames(Retro.protein.all.order)

Retro.new = merge(Retro.protein.all.order, Retro.Figure3.meta.order.2[,c(1,7,9)], by = "Sample")
Retro.new$cluster = factor(Retro.new$hierarchical.clustering, levels = c(3,5,1,2,4))
Retro.new.melt = melt(Retro.new)

theme_set(theme_grey(base_size = 20))
pdf("Retro_Ribosome.pdf", width = 8, height = 4)
ggboxplot(Retro.new.melt[Retro.new.melt$variable %in% ribosome.path,], x = "variable", y = "value", fill = "cluster",
          ylab = "CPTAC relative\nprotein abundance", xlab = "", palette = c("#F8766D","#A3A500","#00BD7A","#00B0F6","#E76BF3"), short.panel.labs = TRUE) +
  theme(text = element_text(size = 16))
dev.off()

theme_set(theme_grey(base_size = 20))
pdf("Retro_EMT.pdf", width = 8, height = 4)
ggboxplot(Retro.new.melt[Retro.new.melt$variable %in% EMT.path,], x = "variable", y = "value", fill = "cluster",
          ylab = "relative protein abundance", xlab = "", palette = c("#F8766D","#A3A500","#00BD7A","#00B0F6","#E76BF3"), short.panel.labs = TRUE) +
  theme(text = element_text(size = 16))
dev.off()

############################################################
## Figure S5A
############################################################
FD.cluster.pca = as.data.frame(t(FD.pathway.order))
FD.cluster.pch = replicate(101, 21)
FD.order.clusters = FD.meta.order$Cluster
FD.order.clusters.1 = rep("white", length(FD.order.clusters))
FD.order.clusters.1[FD.order.clusters=="Cluster_1"] = "#F8766D"
FD.order.clusters.1[FD.order.clusters=="Cluster_2"] = "#A3A500"
FD.order.clusters.1[FD.order.clusters=="Cluster_3"] = "#00BD7A"
FD.order.clusters.1[FD.order.clusters=="Cluster_4"] = "#00B0F6"
FD.order.clusters.1[FD.order.clusters=="Cluster_5"] = "#E76BF3"

FZ.cluster.pca = as.data.frame(t(FZ.pathway.order))
FZ.cluster.pch = replicate(27, 23)

FZ.Figure3.meta.Hierc.NotOverlap.no3 = FZ.Figure3.meta.Hierc.NotOverlap[!FZ.Figure3.meta.Hierc.NotOverlap$hierarchical.clustering == "3",]
FZ.Figure3.meta.Hierc.NotOverlap.no4 = FZ.Figure3.meta.Hierc.NotOverlap.no3[!FZ.Figure3.meta.Hierc.NotOverlap.no3$hierarchical.clustering == "4",]

FZ.cluster.pca.NonOverlap = FZ.cluster.pca[FZ.Figure3.meta.Hierc.NotOverlap.no4$Sample,]
FZ.pathway.order.clusters.Non = c("#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#F8766D","#00B81F","#00B81F","#00B81F","#00B81F","#00B81F","#00B81F","#00B81F","#00B81F","#00B81F","#00B81F","#00B81F","#E76BF3","#E76BF3","#E76BF3","#E76BF3","#E76BF3","#E76BF3")

FD.FZ.cluster.tSNE = rbind(FD.cluster.pca[c(1:33,65:102,129:158),], FZ.cluster.pca.NonOverlap)
FD.FZ.cluster.pch = c(FD.cluster.pch, FZ.cluster.pch)
FD.FZ.cluster.col = c(FD.order.clusters.1[c(1:33,65:102,129:158)], FZ.pathway.order.clusters.Non)

FD.FZ.cluster.tSNE.tSNE = Rtsne(FD.FZ.cluster.tSNE, perplexity = 25, pca_scale = TRUE, check_duplicates = FALSE)

theme_set(theme_grey(base_size = 20))
pdf("tSNE_FD_FZ_Pathway.pdf", height = 8, width = 8)
plot(FD.FZ.cluster.tSNE.tSNE$Y, col = FD.FZ.cluster.col, bg = FD.FZ.cluster.col,
     pch = FD.FZ.cluster.pch, cex = 1.2)
dev.off()

############################################################
## Figure S5D
############################################################
FD.pathway.top.HM.Cluster3 = as.data.frame(FD.pathway.top.HM)
colnames(FD.pathway.top.HM.Cluster3) = c("FD_Cluster1","FD_Cluster2","FD_Cluster3","FD_Cluster4","FD_Cluster5")
FD.pathwayGroups.Cluster3 = FD.Figure3.pathwayGroups
FD.pathwayGroups.Cluster3 = cbind(FD.pathwayGroups.Cluster3, FD.pathway.top.HM.Cluster3)

PDX.pathway = read.table("input_pathway/PDX_pathway_20211216.txt", header = TRUE, quote = "\"", check.names = FALSE, sep = "\t", row.names = 1)
PDX.pathway.Sample10 = PDX.pathway[,grep("Patient_10", colnames(PDX.pathway))]
PDX.pathway.Sample10$Sample10_Average = rowMeans(PDX.pathway.Sample10)
PDX.pathway.Sample10$Pathway = rownames(PDX.pathway.Sample10)
PDX.pathway.Sample10 = merge(PDX.pathway.Sample10, FD.pathwayGroups.Cluster3, by = "Pathway")
PDX.pathway.Sample10$Group = factor(PDX.pathway.Sample10$Group)
PDX.pathway.Sample10$Order = NULL

PDX.pathway.Sample10.melt = melt(PDX.pathway.Sample10, id.vars = c("Pathway","Group"))
PDX.pathway.Sample10.melt = PDX.pathway.Sample10.melt[grep("TMT", PDX.pathway.Sample10.melt$variable),]

PDX.Patient10.melt = ggplot(PDX.pathway.Sample10.melt, aes(x = variable, y = value, col = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(~variable, ncol = 1, scales = "free", drop = TRUE, strip.position = NULL) +
  theme_classic() + labs(x = "", y = "") + theme_bw() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none", text = element_text(size = 14), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

theme_set(theme_grey(base_size = 20))
pdf("PDX_Patient10_Boxplot.pdf", height = 10, width = 4, useDingbats = FALSE)
PDX.Patient10.melt
dev.off()
