############################################################
## title: "Figure S1H_I_J"
## author: "Oscar Murillo"
## date: "August 18, 2020"
############################################################
############################################################
## Figure S1I
############################################################
FD.RNA.quality = read.table("ovarian-dna-qc-formatted.txt", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = NULL)

FD.RNA.quality.meta = merge(FD.protein.quality.meta, FD.RNA.quality, by = "Sample")
FD.RNA.quality.meta$Pool = NULL
FD.RNA.quality.meta$TMT_Plex = NULL
FD.RNA.quality.meta$NonKR = NULL
FD.RNA.quality.meta$Mox = NULL
FD.RNA.quality.meta$MC = NULL
FD.RNA.quality.meta$Cterm = NULL
FD.RNA.quality.meta$RNA_ID = NULL

FD.RNA.quality.meta.melt = melt(FD.RNA.quality.meta)
FD.RNA.quality.meta.melt$Location = gsub("GIOM","Other",FD.RNA.quality.meta.melt$Location)
FD.RNA.quality.meta.melt$Location = gsub("GI","Other",FD.RNA.quality.meta.melt$Location)
FD.RNA.quality.meta.melt$Location = gsub("LN","Other",FD.RNA.quality.meta.melt$Location)
FD.RNA.quality.meta.melt$Location = gsub("Mix","Other",FD.RNA.quality.meta.melt$Location)
FD.RNA.quality.meta.melt$Location = gsub("PT","Other",FD.RNA.quality.meta.melt$Location)
FD.RNA.quality.meta.melt$Location = gsub("UT","Other",FD.RNA.quality.meta.melt$Location)

RNA.quality.keep = c("Number of Reads","Detected Splice Junctions","Duplication Percentage")
FD.RNA.quality.meta.melt = FD.RNA.quality.meta.melt[FD.RNA.quality.meta.melt$variable %in% RNA.quality.keep,]

RNA.response = ggboxplot(FD.RNA.quality.meta.melt, x = "Response", y = "value", ylab = "score", xlab = "tumor response",
                         order = c("Sensitive","Refractory"),
                         palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) + 
  stat_compare_means(comparisons = my_comp_RS, method = "t.test", label = "p.format") +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25)))

RNA.age = ggboxplot(FD.RNA.quality.meta.melt, x = "Age_Category", y = "value", ylab = "score", xlab = "sample age (split by median)",
                    order = c("new","old"),
                    palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) + 
  stat_compare_means(comparisons = my_comp_ON, method = "t.test", label = "p.format") +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25)))

RNA.source = ggboxplot(FD.RNA.quality.meta.melt, x = "Source", y = "value", ylab = "score", xlab = "sample source",
                       order = c("FHCRC","Mayo","UAB"),
                       palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) + 
  stat_compare_means(comparisons = my_comp_UFM, method = "t.test", label = "p.format") +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25)))

theme_set(theme_grey(base_size = 12))
pdf("FD_RNA_Boxplot_All.pdf", width = 8, height = 8)
ggarrange(RNA.response, RNA.age, RNA.source, 
          ncol = 1, nrow = 3)
dev.off()

############################################################
## Figure S1J
############################################################
FD.DNA.quality = read.table("CPTAC3 - Ovaria - DNA WGS Quality.txt", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = NULL)

FD.DNA.quality.meta = merge(FD.DNA.quality, FD.protein.quality.meta, by = "Sample", all.x = TRUE)
FD.DNA.quality.meta$Pool = NULL
FD.DNA.quality.meta$TMT_Plex = NULL
FD.DNA.quality.meta$NonKR = NULL
FD.DNA.quality.meta$Mox = NULL
FD.DNA.quality.meta$MC = NULL
FD.DNA.quality.meta$Cterm = NULL

FD.DNA.quality.meta.melt = melt(FD.DNA.quality.meta)
FD.DNA.quality.meta.melt$Location = gsub("GIOM","Other",FD.DNA.quality.meta.melt$Location)
FD.DNA.quality.meta.melt$Location = gsub("GI","Other",FD.DNA.quality.meta.melt$Location)
FD.DNA.quality.meta.melt$Location = gsub("LN","Other",FD.DNA.quality.meta.melt$Location)
FD.DNA.quality.meta.melt$Location = gsub("Mix","Other",FD.DNA.quality.meta.melt$Location)
FD.DNA.quality.meta.melt$Location = gsub("PT","Other",FD.DNA.quality.meta.melt$Location)
FD.DNA.quality.meta.melt$Location = gsub("UT","Other",FD.DNA.quality.meta.melt$Location)

DNA.quality.keep = c("TOTAL_DUPLICATION","PCT_CHIMERAS","GC_DROPOUT")
FD.DNA.quality.meta.melt = FD.DNA.quality.meta.melt[FD.DNA.quality.meta.melt$variable %in% DNA.quality.keep,]

DNA.response = ggboxplot(FD.DNA.quality.meta.melt, x = "Response", y = "value", ylab = "score", xlab = "tumor response",
                         order = c("Sensitive","Refractory"),
                         palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) + 
  stat_compare_means(comparisons = my_comp_RS, method = "t.test", label = "p.format") +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25)))

DNA.age = ggboxplot(FD.DNA.quality.meta.melt, x = "Age_Category", y = "value", ylab = "score", xlab = "sample age (split by median)",
                    order = c("new","old"),
                    palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) + 
  stat_compare_means(comparisons = my_comp_ON, method = "t.test", label = "p.format") +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25)))

DNA.source = ggboxplot(FD.DNA.quality.meta.melt, x = "Source", y = "value", ylab = "score", xlab = "sample source",
                       order = c("FHCRC","Mayo","UAB"),
                       palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) + 
  stat_compare_means(comparisons = my_comp_UFM, method = "t.test", label = "p.format") +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25)))

theme_set(theme_grey(base_size = 12))
pdf("FD_DNA_Boxplot_All.pdf", width = 8, height = 8)
ggarrange(DNA.response, DNA.age, DNA.source, 
          ncol = 1, nrow = 3)
dev.off()

############################################################
## Figure S1H
############################################################
FD.protein.quality = read.table("quality/Protein_Quality_20211018.txt", header = TRUE, quote="\"", check.names = FALSE, sep = "\t", row.names = NULL)
FD.protein.quality = FD.protein.quality[,1:12]

FD.protein.quality.meta = merge(FD.meta, FD.protein.quality, by = "Sample")
FD.protein.quality.meta$Age = NULL
FD.protein.quality.meta$Grade = NULL
FD.protein.quality.meta$Age_Category = ifelse(FD.protein.quality.meta$Sample_Age > median(FD.protein.quality.meta$Sample_Age), "old","new")
FD.protein.quality.meta$TMT_ratio = NULL
FD.protein.quality.meta$ALB_ratio = NULL
FD.protein.quality.meta$HB_ratio = NULL
FD.protein.quality.meta$Hist_ratio = NULL
FD.protein.quality.meta$Conc = NULL
FD.protein.quality.meta$Sample_Age = NULL
FD.protein.quality.melt.2 = melt(FD.protein.quality.meta)

FD.protein.quality.melt.2$Location = gsub("GIOM","Other",FD.protein.quality.melt.2$Location)
FD.protein.quality.melt.2$Location = gsub("GI","Other",FD.protein.quality.melt.2$Location)
FD.protein.quality.melt.2$Location = gsub("LN","Other",FD.protein.quality.melt.2$Location)
FD.protein.quality.melt.2$Location = gsub("Mix","Other",FD.protein.quality.melt.2$Location)
FD.protein.quality.melt.2$Location = gsub("PT","Other",FD.protein.quality.melt.2$Location)
FD.protein.quality.melt.2$Location = gsub("UT","Other",FD.protein.quality.melt.2$Location)

FD.protein.quality.melt.2$variable = gsub("NonKR","SemiKR",FD.protein.quality.melt.2$variable)
FD.protein.quality.melt.2$variable = factor(FD.protein.quality.melt.2$variable, levels = c("SemiKR","Cterm","MC","Mox"))

Protein.Response = ggboxplot(FD.protein.quality.melt.2, x = "Response", y = "value", ylab = "score", xlab = "tumor response",
                             palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) + 
  stat_compare_means(comparisons = my_comp_RS, method = "t.test", label = "p.format", label.y = c(1.5)) +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.7,1.75), expand = expansion(mult = c(0.1, 0.1)))

Protein.age = ggboxplot(FD.protein.quality.melt.2, x = "Age_Category", y = "value", ylab = "score", xlab = "sample age (split by median)",
                        order = c("new","old"),
                        palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) + 
  stat_compare_means(comparisons = my_comp_ON, method = "t.test", label = "p.format", label.y = c(1.5)) +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.7,1.75), expand = expansion(mult = c(0.1, 0.1)))

Protein.source = ggboxplot(FD.protein.quality.melt.2, x = "Source", y = "value", ylab = "score", xlab = " sample source",
                           palette = "jco", add = "jitter", facet.by = "variable", short.panel.labs = FALSE) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) + 
  stat_compare_means(comparisons = my_comp_UFM, method = "t.test", label = "p.format", label.y = c(1.5, 1.6, 1.7)) +
  theme(text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.7,1.75), expand = expansion(mult = c(0.1, 0.1)))

theme_set(theme_grey(base_size = 12))
pdf("FD_Protein_Boxplot_All.pdf", width = 8, height = 8)
ggarrange(Protein.Response, Protein.age, Protein.source, 
          ncol = 1, nrow = 3)
dev.off()