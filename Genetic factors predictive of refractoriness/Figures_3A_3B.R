## Figure 3-A
```{r}
library(rstatix)
library(data.table)
library(BSgenome)

https://urldefense.proofpoint.com/v2/url?u=http-3A__seg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=8KkBanu2jNwwU_b2gdyK1e88IM_p29x2EQ6kTrL9O6A&e=  <- fread("CNV segmentation data")
meta             <- fread("Meta data")
cnvex.purity     <- fread("Purity ploidy measures file") ## we want purity values from this file
https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e=        <- fread("CHR17 LOH data") ## using this list of samples
hg38.cyto.coords  <- fread("HG38 cytobands downloaded as GRanges")
colnames(cnvex.purity) <- c("sample", "cand", "CNVEX.ploidy", "CNVEX.purity")

meta <- meta[patient %in% chr17.info$case]
meta <- merge.data.table(meta, cnvex.purity, all.x = T, by.x = "DNA ID", by.y = "sample")


## order arms
tmp <- findOverlaps(makeGRangesFromDataFrame(https://urldefense.proofpoint.com/v2/url?u=http-3A__seg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=8KkBanu2jNwwU_b2gdyK1e88IM_p29x2EQ6kTrL9O6A&e= ), hg38.cyto.coords, select="first")
seg.data$cytoband <- hg38.cyto.coords[tmp]$name
seg.data$arm <- paste0(seg.data$seqnames, str_sub(seg.data$cytoband, 1, 1))
seg.data$arm <- factor(seg.data$arm, unique(seg.data$arm), ordered=TRUE)

#### start the heatmap
oc.age = meta[, .(id=patient, age=PatientAge)]
oc.age = meta[, .(id=patient, NeoAdjuvant=NeoAdjuvant)]

fns <- data.table(digest.fn = digest.fns)
fns[, ids := str_replace(basename(digest.fn), "-somatic-digest.rds", "")]
fns[, patient := tstrsplit(ids,"-")[1]]
fns <- fns[patient %in% chr17.info$case]
doParallel::registerDoParallel(30)
arm.level.prev <- foreach(i = 1:nrow(meta)) %dopar% {
  id <- meta[i,]$patient
  seg <- seg.data[dna_id == id]
  
  meta.sample <- meta[i, ]
  
  arm.level.prev <- seg[, .(prev = sum((C != weighted.median(C, width) | K != weighted.median(K, width)) & C != meta.sample$median.C, na.rm = TRUE)/(length(C))), by = .(arm)] ## all # not equal to sample, because HRD might cause >50% change
  arm.level.prev$case = id
  ## arm level
  return(as.data.table(arm.level.prev))
}

arm.level.prev <- rbindlist(arm.level.prev)
focal.prevelance <- arm.level.prev[, .(perc = mean(prev, na.rm=TRUE)),  by = case]

focal.prevelance <- merge.data.table(focal.prevelance, meta[, c("patient", "Response")], by.y = "patient", by.x = "case")


stat.focal <- as.data.frame(focal.prevelance) %>%
  wilcox_test(perc ~ Response) %>%
  add_significance() %>% add_xy_position()
stat.focal$y.position <- 0.33
p.focal = focal.prevelance %>% ggplot(aes(x = Response, y = perc))+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_quasirandom(color = "#1ce0b2", varwidth = TRUE, width = 0.21, size = 4)+
  stat_pvalue_manual(stat.focal, label = "P = {p}",vjust = 0.1, bracket.nudge.y = -0.01, size = 7)+
  theme_pubr()+
  labs(y = "% Genome w/ sub-arm events")+
  theme(axis.title = element_text(size = 25), axis.text =  element_text(size = 18))+
  scale_y_continuous(, breaks = c(0, 0.1, 0.2, 0.3), labels = scales::percent(c(0, 0.1, 0.2, 0.3), accuracy = 1))

p.focal
```

## Figure 3-B
```{r}
tp53_scores      <- fread("TP53 scores")
meta             <- fread("Meta data")
cnvex.purity     <- fread("Purity ploidy measures file") ## we want purity values from this file
https://urldefense.proofpoint.com/v2/url?u=http-3A__cin.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=mqYaOqowyNIp7rgA9YHdivhfaBUOt4uFykZc7QJvnDE&e=          <- fread("CIN scores")
mutations.strict <- fread("Mutation data")
https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e=        <- fread("CHR17 LOH data") ## using this list of samples
mutation.tbl     <- fread("FFPE mutation")
tmb              <- fread("TMB data")
colnames(meta)[1] = "patient"
colnames(cnvex.purity) <- c("sample", "cand", "CNVEX.ploidy", "CNVEX.purity")

meta <- meta[patient %in% chr17.info$case]
meta <- merge.data.table(meta, cnvex.purity, all.x = T, by.x = "patient", by.y = "sample")

https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.group&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=IsesMnh5lXnj1RF1MpLfZbvqd4xH5JXITeb71RDjIWc&e=  = chr17.info[`Chr17-LOH` == 1]$case
https://urldefense.proofpoint.com/v2/url?u=http-3A__other.group&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=oAKXwLWrscIHdykxoAvua9MuDotoPElOHpxXGdOw3cQ&e=  = chr17.info[`Chr17-LOH` == 0]$case

## fix mutation table colnames
mut.names <- colnames(mutation.tbl)[-1]
name.tbl  <- meta[, c("PatientID", "patient")]
setkey(name.tbl, PatientID)
name.tbl  <- name.tbl[mut.names]
colnames(mutation.tbl) <- c("Gene_ID", name.tbl$patient)

## mutations pick
tp53.mutation <- t(mutation.tbl[Gene_ID == "TP53"][,-1])
tp53.mutation <- data.table(sample = rownames(tp53.mutation), mutation = as.vector(tp53.mutation))

brca.mutation <- t(mutation.tbl[Gene_ID %in% c("BRCA1", "BRCA2")][,-1])
brca.mutation.comb <- as.numeric(brca.mutation[,1] + brca.mutation[,2] > 0)
brca.mutation <- data.table(sample = rownames(brca.mutation), mutation = as.vector(brca.mutation.comb))

doParallel::registerDoParallel(30)

https://urldefense.proofpoint.com/v2/url?u=http-3A__loh.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=g_yZYstD33Zi72tHhG0U60un8tVrTvnT0gXr8srKdhw&e=  <- foreach(i = 1:nrow(meta)) %dopar% {
  seg <-  seg.data[dna_id == meta[i,]$patient]
  
  loh       <- weighted.mean(seg$K == 0 | seg$C == 1, seg$width, na.rm=TRUE)
  return(data.table(loh = loh, case = meta[i,]$patient))
}
https://urldefense.proofpoint.com/v2/url?u=http-3A__loh.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=g_yZYstD33Zi72tHhG0U60un8tVrTvnT0gXr8srKdhw&e=  <- rbindlist(https://urldefense.proofpoint.com/v2/url?u=http-3A__loh.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=g_yZYstD33Zi72tHhG0U60un8tVrTvnT0gXr8srKdhw&e= )

https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- data.table(case = https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.group&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=IsesMnh5lXnj1RF1MpLfZbvqd4xH5JXITeb71RDjIWc&e= , group= "noLOH")
https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- rbind(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , data.table(case = https://urldefense.proofpoint.com/v2/url?u=http-3A__other.group&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=oAKXwLWrscIHdykxoAvua9MuDotoPElOHpxXGdOw3cQ&e= , group= "LOH"))
https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , meta[, c("patient", "tAi5", "Response")], by.x = "case", by.y = "patient")
https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , https://urldefense.proofpoint.com/v2/url?u=http-3A__loh.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=g_yZYstD33Zi72tHhG0U60un8tVrTvnT0gXr8srKdhw&e= , by = "case")
https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , https://urldefense.proofpoint.com/v2/url?u=http-3A__cin.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=mqYaOqowyNIp7rgA9YHdivhfaBUOt4uFykZc7QJvnDE&e= , by.x = "case", by.y = "sample")

reg.data$tp53_mut <- as.numeric(reg.data$case %in% tp53.mutation[mutation == 1]$sample)
reg.data$brca_mut <- as.numeric(reg.data$case %in% brca.mutation[mutation == 1]$sample)

https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- reg.data[case %in% chr17.info$case]
reg.data$chr17LOH <- ifelse(reg.data$group == "LOH", 1 , 0)
reg.data$tAi5 = as.numeric(reg.data$tAi5)
reg.data$loh.high <- ifelse(reg.data$loh > 0.15, 1, 0)
reg.data$Response <- as.factor(reg.data$Response)

tmb <- merge.data.table(tmb, meta[, c("PatientID", "patient")], by.x = "patient_id", by.y = "PatientID")
colnames(tmb)[2] <- "Mutation count"
colnames(tmb)[4] <- "case"
https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , tmb[, c("case", "Mutation count")])

library(forestmodel)
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "tAi5"] <- "nTAI"
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "loh"] <- "%LOH"
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "chr17LOH"] <- "Chr17-LOH"
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "loh15"] <- "nLOH"
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "lst10"] <- "nLST"
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "ploidy"] <- "Ploidy"
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "tp53_mut"] <- "TP53 Mutation"
colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )[colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= ) == "brca_mut"] <- "BRCA1/2 Mutation"

reg.data$nTAI <- scale(reg.data$nTAI)
reg.data$nLST <- scale(reg.data$nLST)
reg.data$nLOH <- scale(reg.data$nLOH)
reg.data$`Mutation count` <- scale(reg.data$`Mutation count`)
reg.data$Ploidy <- scale(reg.data$Ploidy)

fit <- glm(Response ~ `Chr17-LOH`+nTAI+nLST+nLOH+`TP53 Mutation`+`BRCA1/2 Mutation`+`Mutation count`, data = https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , family = binomial)

summary(fit)

format = forest_model_format_options()
format$text_size = 8
format$point_size = 12
format$shape = 12
format$colour = "red"
  forest_model(fit, format_options = format)
  
  ```
  
  
  