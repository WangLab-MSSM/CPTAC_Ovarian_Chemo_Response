## Figure 2-C:
```{r}
library(rstatix)
library(data.table)

https://urldefense.proofpoint.com/v2/url?u=http-3A__seg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=8KkBanu2jNwwU_b2gdyK1e88IM_p29x2EQ6kTrL9O6A&e=  <- fread("CNV segmentation data")
meta             <- fread("Meta data")
cnvex.purity     <- fread("Purity ploidy measures file") ## we want purity values from this file
https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e=        <- fread("CHR17 LOH data") ## using this list of samples

colnames(cnvex.purity) <- c("sample", "cand", "CNVEX.ploidy", "CNVEX.purity")

meta <- meta[patient %in% chr17.info$case]
meta <- merge.data.table(meta, cnvex.purity, all.x = T, by.x = "DNA ID", by.y = "sample")



thr = 0.25
https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e=  <- foreach(i = 1:nrow(https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e= )) %dopar% {
  id <- chr17.info[i, ]$case
  seg       <- seg.data[dna_id == id]
  median.C  <- weighted.median(as.numeric(seg$C), as.numeric(seg$width), na.rm = TRUE)
  C.dosage  <- seg$C - median.C
  K         <- seg$K
  lr        <- seg$lr
  return(data.table(C.dosage = C.dosage, K = K, sample = id, lr = lr, width = seg$width,
                    chr = as.character(seg$seqnames), idx = 1:nrow(seg)))
}
https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e=  <- rbindlist(https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e= )
chr.data.loh <- cnv.data[ , .(loh  = weighted.mean(K == 0, width, na.rm=TRUE)), by = .(sample, chr)]
chr.data.gai <- cnv.data[ , .(gain = weighted.mean(C.dosage > 0, width, na.rm=TRUE)), by = .(sample, chr)]
chr.data.los <- cnv.data[ , .(loss = weighted.mean(C.dosage < 0, width, na.rm=TRUE)), by = .(sample, chr)]
https://urldefense.proofpoint.com/v2/url?u=http-3A__chr.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=02A0cZVCZoHVkRYuzem1Chf7bPY6jHNyzLpE5GnP854&e=      <- cbind(chr.data.loh,chr.data.gai[, "gain"], chr.data.los[, "loss"])

https://urldefense.proofpoint.com/v2/url?u=http-3A__chr.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=02A0cZVCZoHVkRYuzem1Chf7bPY6jHNyzLpE5GnP854&e=  <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__chr.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=02A0cZVCZoHVkRYuzem1Chf7bPY6jHNyzLpE5GnP854&e= , meta[, c("DNA ID", "Response", "tAi5", "CNVEX.purity")], all.x=TRUE, by.x = "sample", by.y = "DNA ID")

## all the other loh's
loh.pvals <- lapply(paste0("chr", 1:22), function(select.chr) {
  
  https://urldefense.proofpoint.com/v2/url?u=http-3A__onlychr.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=A-byACUhSgJIWD2AUk8M9PsCJWjrnzLdn_SFCyZ1oGQ&e=  <- list("low"  = split(chr.data[chr %in% select.chr], chr.data[chr %in% select.chr]$loh <= thr)[[2]],
                                                                                                                                                                                                                                                                                                     "high" = split(chr.data[chr %in% select.chr], chr.data[chr %in% select.chr]$loh >= 0.5)[[2]])
  https://urldefense.proofpoint.com/v2/url?u=http-3A__test.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=uyJ6s-m13NPOU-_EAmOfr2zEAXE-H9sLlgyw38e3x7g&e=       <- matrix(c(sum(onlychr.data[[1]]$Response == "Sensitive"), sum(onlychr.data[[2]]$Response == "Sensitive"),
                                                                                                                                                                                                                                                                                                         sum(onlychr.data[[1]]$Response == "Refractory"), sum(onlychr.data[[2]]$Response == "Refractory")), nrow=2)
  fisher.p     <- fisher.test(https://urldefense.proofpoint.com/v2/url?u=http-3A__test.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=uyJ6s-m13NPOU-_EAmOfr2zEAXE-H9sLlgyw38e3x7g&e= )$p.value
  return(data.table(type = "loh",chr = select.chr, fisher.p = fisher.p))
})

loh.pvals <- rbindlist(loh.pvals)
loh.pvals <- adjust_pvalue(loh.pvals, "fisher.p")
## all the other gains
gain.pvals <- lapply(paste0("chr", 1:22), function(select.chr) {
  # https://urldefense.proofpoint.com/v2/url?u=http-3A__onlychr.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=A-byACUhSgJIWD2AUk8M9PsCJWjrnzLdn_SFCyZ1oGQ&e=  <- split(chr.data[chr %in% select.chr], chr.data[chr %in% select.chr]$gain <= thr)
  https://urldefense.proofpoint.com/v2/url?u=http-3A__onlychr.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=A-byACUhSgJIWD2AUk8M9PsCJWjrnzLdn_SFCyZ1oGQ&e=  <- list("low"  = split(chr.data[chr %in% select.chr], chr.data[chr %in% select.chr]$gain <= thr)[[2]], ## less than 25%
                                                                                                                                                                                                                                                                                                     "high" = split(chr.data[chr %in% select.chr], chr.data[chr %in% select.chr]$gain >= 0.5)[[2]]) ## more than half
  https://urldefense.proofpoint.com/v2/url?u=http-3A__test.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=uyJ6s-m13NPOU-_EAmOfr2zEAXE-H9sLlgyw38e3x7g&e=       <- matrix(c(sum(onlychr.data[[1]]$Response == "Sensitive"), sum(onlychr.data[[2]]$Response == "Sensitive"),
                                                                                                                                                                                                                                                                                                         sum(onlychr.data[[1]]$Response == "Refractory"), sum(onlychr.data[[2]]$Response == "Refractory")), nrow=2)
  fisher.p     <- fisher.test(https://urldefense.proofpoint.com/v2/url?u=http-3A__test.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=uyJ6s-m13NPOU-_EAmOfr2zEAXE-H9sLlgyw38e3x7g&e= )$p.value
  return(data.table(type = "gain",chr = select.chr, fisher.p = fisher.p))
})
gain.pvals <- rbindlist(gain.pvals)
gain.pvals <- adjust_pvalue(gain.pvals, "fisher.p")
## all the other deletions
loss.pvals <- lapply(paste0("chr", 1:22), function(select.chr) {
  https://urldefense.proofpoint.com/v2/url?u=http-3A__onlychr.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=A-byACUhSgJIWD2AUk8M9PsCJWjrnzLdn_SFCyZ1oGQ&e=  <- list("low"  = split(chr.data[chr %in% select.chr], chr.data[chr %in% select.chr]$loss <= thr)[[2]],
                                                                                                                                                                                                                                                                                                     "high" = split(chr.data[chr %in% select.chr], chr.data[chr %in% select.chr]$loss >= 0.5)[[2]])
  https://urldefense.proofpoint.com/v2/url?u=http-3A__test.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=uyJ6s-m13NPOU-_EAmOfr2zEAXE-H9sLlgyw38e3x7g&e=       <- matrix(c(sum(onlychr.data[[1]]$Response == "Sensitive"), sum(onlychr.data[[2]]$Response == "Sensitive"),
                                                                                                                                                                                                                                                                                                         sum(onlychr.data[[1]]$Response == "Refractory"), sum(onlychr.data[[2]]$Response == "Refractory")), nrow=2)
  fisher.p     <- fisher.test(https://urldefense.proofpoint.com/v2/url?u=http-3A__test.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=uyJ6s-m13NPOU-_EAmOfr2zEAXE-H9sLlgyw38e3x7g&e= )$p.value
  return(data.table(type = "loss",chr = select.chr, fisher.p = fisher.p))
})
loss.pvals <- rbindlist(loss.pvals)
loss.pvals <- adjust_pvalue(loss.pvals, "fisher.p")

pvals <- rbind(loh.pvals, gain.pvals, loss.pvals)
pvals$chr <- factor(pvals$chr, levels = unique(pvals$chr))
## plot

pvals2 <- copy(pvals)
pvals2$type <- fcase(pvals2$type == "loss", "LOSS",
                     pvals2$type == "loh" , "LOH" ,
                     pvals2$type == "gain", "GAIN")
pvals2.base = data.table(type = rep("BASE",table(pvals2$type)[1]), chr = pvals2[type == "LOH"]$chr, fisher.p = 1, fisher.p.adj = 1)
pvals2 <- rbind(pvals2, pvals2.base)

pvals2[, col := fcase(type == "LOH", "#69ef7b",
                      type == "GAIN", "#923124",
                      type == "LOSS", "#79e3f9",
                      type == "BASE", "white")]
pvals2[, size := fcase(type == "LOH", 2,
                       type == "GAIN", 2,
                       type == "LOSS", 2,
                       type == "BASE", 0)]
pvals2[, alpha := fcase(type == "LOH", 1,
                        type == "GAIN", 1,
                        type == "LOSS", 1,
                        type == "BASE", 0)]
pvals2[, shape := fcase(type == "LOH", 15,
                        type == "GAIN",17,
                        type == "LOSS",16)]


pvals2$chr <- factor(str_replace(pvals2$chr, "chr", ""), levels = 1:22)
p1 = pvals2 %>% ggplot()+
  geom_line(aes(x = chr,y = -log10(fisher.p)))+
  geom_point(aes(x = chr, y = -log10(fisher.p), alpha = alpha), size = 5, color = pvals2$col)+
  theme_pubr()+labs(x = "Chromosome")+geom_hline(yintercept = -log10(0.05), linetype = "dotted")+
  theme(axis.title = element_blank(), axis.text =  element_text(size = 18), legend.position = "none")
p2 = pvals2 %>% ggplot()+geom_line(aes(x = chr, y = -log10(fisher.p.adj)))+
  geom_point(aes(x = chr, y = -log10(fisher.p.adj), alpha = alpha), size = 5, color = pvals2$col, shape = pvals2$shape)+
  labs(x = "Chromosome", y = "-Log10(q-value)")+
  theme_pubr()+geom_hline(yintercept = -log10(0.1), linetype = "dotted")+theme(legend.position = "none")+
  theme(axis.title = element_blank(), axis.text =  element_text(size = 18), legend.position = "none")

# p1 / p2
lgd.qvalues = ComplexHeatmap::Legend(at = c("LOSS","GAIN","LOH"),
                                     labels = c("Deletion", "Amplification", "Loss of heterozygousity"),
                                     legend_gp = gpar(fill=c("LOSS" = "#79e3f9","GAIN" = "#923124","LOH" = "#69ef7b")),
                                     title = "Event",
                                     nrow = 1,
                                     graphics = list(
                                       function(x, y, w, h) grid.points(x, y, gp = gpar(col = "#79e3f9"), pch = 16),
                                       function(x, y, w, h) grid.points(x, y, gp = gpar(col = "#923124"), pch = 17),
                                       function(x, y, w, h) grid.points(x, y, gp = gpar(col = "#69ef7b"), pch = 15)
                                     )
)


```

## Figure 2-E:
```{r}
library(survival)
library(survminer)
library(data.table)

https://urldefense.proofpoint.com/v2/url?u=http-3A__sample.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=xQtiaPw3qaK5GvwFL1-veouqDnrFZJFYptTBoVfkxNg&e=    <- fread("MSK IMPACT clinical sample data (data_clinical_sample.txt)")
https://urldefense.proofpoint.com/v2/url?u=http-3A__patient.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OWrtfgIGpHXCcmx5NMCL8kIZdID6TQR1hEEKEJUGgDU&e=   <- fread("MSK IMPACT clinical patient data (data_clinical_patient.txt)")
https://urldefense.proofpoint.com/v2/url?u=http-3A__cna.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=GpOr3UwKTXbFx4mxyYpVHHgs7-GPbM0C-qRvnmhXKiw&e=       <- fread("MSK IMPACT CNA data (data_cna.txt)")

https://urldefense.proofpoint.com/v2/url?u=http-3A__sample.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=xQtiaPw3qaK5GvwFL1-veouqDnrFZJFYptTBoVfkxNg&e=    <- sample.data[`Subtype Abbreviation` == "HGSOC",] 


https://urldefense.proofpoint.com/v2/url?u=http-3A__patient.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OWrtfgIGpHXCcmx5NMCL8kIZdID6TQR1hEEKEJUGgDU&e=   <- patient.data[`#Patient Identifier` %in% sample.data$`Patient Identifier`]
https://urldefense.proofpoint.com/v2/url?u=http-3A__patient.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OWrtfgIGpHXCcmx5NMCL8kIZdID6TQR1hEEKEJUGgDU&e=   <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__patient.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OWrtfgIGpHXCcmx5NMCL8kIZdID6TQR1hEEKEJUGgDU&e= , sample.data[, c("Patient Identifier", "#Sample Identifier")], 
                                                                                                                                                                                                                                                                                                                by.x = "#Patient Identifier", by.y = "Patient Identifier")

samples.col   <- which(colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__cna.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=GpOr3UwKTXbFx4mxyYpVHHgs7-GPbM0C-qRvnmhXKiw&e= ) %in% sample.data$`#Sample Identifier`)
https://urldefense.proofpoint.com/v2/url?u=http-3A__cna.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=GpOr3UwKTXbFx4mxyYpVHHgs7-GPbM0C-qRvnmhXKiw&e=       <- cna.data[,.SD,.SDcols=c(1,samples.col)]
cna.data.tp53 <- cna.data[Hugo_Symbol == "TP53"]


## looking at 17 loss
https://urldefense.proofpoint.com/v2/url?u=http-3A__cna.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=GpOr3UwKTXbFx4mxyYpVHHgs7-GPbM0C-qRvnmhXKiw&e=    <- cna.data[ID %in% sample.data$`#Sample Identifier` & chrom == 17,] 
https://urldefense.proofpoint.com/v2/url?u=http-3A__cna.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=GpOr3UwKTXbFx4mxyYpVHHgs7-GPbM0C-qRvnmhXKiw&e=    <- cna.data[, .(https://urldefense.proofpoint.com/v2/url?u=http-3A__mean.lr&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=SZD9uPIFH_BYk7cGJNJ3XRV-UQ-RoKZxwu7IuVZRaCU&e=  = weighted.mean(seg.mean, loc.end-loc.start)), by = .(ID)]

ggplot()+geom_density(aes(x = cna.data$https://urldefense.proofpoint.com/v2/url?u=http-3A__mean.lr&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=SZD9uPIFH_BYk7cGJNJ3XRV-UQ-RoKZxwu7IuVZRaCU&e= ))+theme_pubr()+scale_x_continuous(breaks = seq(-1,1, 0.1))+
  geom_vline(xintercept = -0.1138661, color = "red")+
  geom_vline(xintercept = 0.04613389, color = "red")

## define chr17 del,wt samples according to the thr
del.samples <- patient.data[`#Sample Identifier` %in% cna.data[https://urldefense.proofpoint.com/v2/url?u=http-3A__mean.lr&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=SZD9uPIFH_BYk7cGJNJ3XRV-UQ-RoKZxwu7IuVZRaCU&e=  < -0.1138661]$ID]
wt.samples  <- patient.data[(`#Sample Identifier` %in% cna.data[https://urldefense.proofpoint.com/v2/url?u=http-3A__mean.lr&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=SZD9uPIFH_BYk7cGJNJ3XRV-UQ-RoKZxwu7IuVZRaCU&e=  > 0.04613389]$ID)]

chr17.status      <- data.table(group = c(rep("chr17-Gain", nrow(wt.samples)), rep("chr17-Loss", nrow(del.samples))),
                                os    = c(wt.samples$`Overall Survival (Months)`, del.samples$`Overall Survival (Months)`),
                                os.s  = c(wt.samples$`Overall Survival Status`, del.samples$`Overall Survival Status`),
                                sample.indentifier = c(wt.samples$`#Sample Identifier`, del.samples$`#Sample Identifier`))

chr17.status$os.s <- as.numeric(tstrsplit(tp53$os.s, ":")[[1]])
chr17.status$os   <- as.numeric(tp53$os)

chr17.status$group <- ifelse(chr17.status$group == "chr17-Gain", "17-Het", "17-LOH")
fit <- survfit(Surv(os, os.s) ~ group, data = chr17.status)

p = ggsurvplot(fit,
               pval = TRUE, https://urldefense.proofpoint.com/v2/url?u=http-3A__conf.int&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=UXvLDpy_Mx_jesrUiaT6_wE8tBAW4b4OJwJAKm8T6-U&e=  = TRUE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               # palette = c("#E7B800", "#2E9FDF")
               font.title = c(26, "black"),
               font.subtitle = c(25, "black"),
               font.caption = c(24, "plain", "black"),
               font.x = c(24, "plain", "black"),
               font.y = c(24, "plain", "black"),
               font.tickslab = c(22, "plain", "black")
)
p$table <- p$table+theme(text = element_text(size = 22))
p$table$layers[[1]]$aes_params$size = 7
p


```

## Figure S2B:
```{r}
## this step needs raw data downloaded from dbgap and ran CNVEX using default setting
library(gtrellis)
rna.oc       <- fread("Normalized RNA expression data")
prot.oc      <- fread("Normalized RNA expression data")
prot.oc      <- fread("Chr17 LOH data")
meta         <- fread("Meta data")
https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e=      <- fread("Gene level Copy number values with median CN per sample moved to 0")
cnvex.purity <- fread("Purity ploidy measures file") ## we want purity values from this file

colnames(cnvex.purity) <- c("sample", "cand", "CNVEX.ploidy", "CNVEX.purity")

meta <- meta[patient %in% chr17.info$case]
meta <- merge.data.table(meta, cnvex.purity, all.x = T, by.x = "DNA ID", by.y = "sample")


## removing the pilot samples:

cohort <- readRDS("Reading all cnvex digest files")
cohort <- merge.data.table(cohort, meta[, c("DNA ID", "Cluster")], by = "DNA ID")
cohort <- merge.data.table(cohort, meta[, c("DNA ID", "Response")], by = "DNA ID")
cohort.split <- split(cohort, cohort$Cluster)

plotCohortDosageFrequency <- function(cohort) {
  chorot.size <- length(cohort)
  https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e=  <- foreach(i = 1:length(cohort)) %dopar% {
    digest    <- readRDS(cohort[i])
    seg       <- cnvex::segOut(digest)
    median.C  <- weighted.median(as.numeric(seg$C), as.numeric(width(seg)), na.rm = TRUE)
    C.dosage  <- digest$tile$C - median.C
    K         <- digest$tile$K   
    return(data.table(C.dosage = C.dosage, K = K, sample = str_replace(basename(cohort), "-wgs-somatic-digest.rds", "")[i], 
                      chr = as.character(seqnames(digest$tile)), idx = 1:length(digest$tile)))
  }
  https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e=  <- rbindlist(https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e= )
  
  C.dosage <- dcast.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e= ,  idx ~ sample, value.var = "C.dosage")
  C.dosage <- C.dosage[, -1]
  gains <- lapply(1:ncol(C.dosage), function(i) {
    # col <- C.dosage[,i, with=FALSE]
    col   <- C.dosage[,get(names(C.dosage)[i])]
    col[is.na(col)]   <- 0
    col.g <- ifelse(col > 0, 1, 0)
    return(col.g)
  })
  gains <- do.call(cbind, gains)
  
  loss <- lapply(1:ncol(C.dosage), function(i) {
    # col <- C.dosage[,i, with=FALSE]
    col   <- C.dosage[,get(names(C.dosage)[i])]
    col[is.na(col)]   <- 0
    col.g <- ifelse(col < 0, 1, 0)
    return(col.g)
  })
  loss <- do.call(cbind, loss)
  
  K.tbl <- dcast.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e= ,  idx ~ sample, value.var = "K")
  K.tbl <- K.tbl[, -1]
  loh <- lapply(1:ncol(K.tbl), function(i) {
    col   <- K.tbl[,get(names(K.tbl)[i])]
    col[is.na(col)]   <- 0
    col.g <- ifelse(col == 0, 1, 0)
    return(col.g)
  })
  loh <- do.call(cbind, loh)
  
  digest    <- readRDS(cohort[1])
  return(data.table(chr = as.character(seqnames(digest$tile)), start= as.numeric(start(digest$tile)), end= as.numeric(end(digest$tile)),
                    loss = rowMeans(loss), gain = rowMeans(gains), loh = rowMeans(loh)))
}
doParallel::registerDoParallel(30)

https://urldefense.proofpoint.com/v2/url?u=http-3A__plot.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=lpmekWdTh4uiS0siIPcWrYMWzvT9WWdF1Vqq1uQUIkY&e=  <- plotCohortDosageFrequency(cohort$digest.fn)
plot.data$loh <- ifelse(plot.data$loh > 0.99, 0, plot.data$loh)

## gtrellis plot!
https://urldefense.proofpoint.com/v2/url?u=http-3A__plot.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=lpmekWdTh4uiS0siIPcWrYMWzvT9WWdF1Vqq1uQUIkY&e=  <- plot.data[chr %in% paste0("chr", 1:22)]
plot.data$chr <- str_replace(plot.data$chr, "chr", "")

gtrellis_layout(byrow = FALSE, n_track = 2, ncol = length(unique(plot.data$chr)),
                add_ideogram_track = TRUE, remove_chr_prefix = TRUE, 
                track_ylim = c( 0 , 1, 
                                0 , 1),
                track_axis = c(FALSE,TRUE),
                track_height = unit.c(
                  2*grobHeight(textGrob("22")),
                  unit(5, "null")),
                category = paste0("chr", 1:22),
                track_ylab = c("","% LOH"),
                border = FALSE)

add_track(panel_fun = function(gr) {
  # the use of `get_cell_meta_data()` will be introduced later
  chr = get_cell_meta_data("name")  
  grid.rect(gp = gpar(fill = "#EEEEEE"))
  grid.text(chr, gp = gpar(fontsize = 7))
})


add_lines_track(as.data.frame(plot.data[,1:3]), plot.data[[6]], area = TRUE, baseline = 0,
                gp = gpar(fill = "#add465"))


```

