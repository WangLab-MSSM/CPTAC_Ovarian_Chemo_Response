## Figure S2C:
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
chr.data.loh <- cnv.data[ , .(loh  = weighted.mean(K == 0, width,na.rm=TRUE)), by = .(sample, chr)]
chr.data.gai <- cnv.data[ , .(gain = weighted.mean(C.dosage > 0, width, na.rm=TRUE)), by = .(sample, chr)]
chr.data.los <- cnv.data[ , .(loss = weighted.mean(C.dosage < 0, width, na.rm=TRUE)), by = .(sample, chr)]
https://urldefense.proofpoint.com/v2/url?u=http-3A__chr.data.lr&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=_N1pMA5M80Z5pmj-OprJawNNgxjrV5tmehnqC5A1BlQ&e=   <- cnv.data[ , .(lr = weighted.mean(lr, width, na.rm=TRUE)), by = .(sample, chr)]


chr.data.lr.loh <- chr.data.lr[chr == "chr17" & sample %in% chr17.info[group == "LOH"]$case]
chr.data.lr.loh$type = "LOH"
chr.data.lr.het <- chr.data.lr[chr == "chr17" & sample %in% chr17.info[group != "LOH"]$case]
chr.data.lr.het$type = "Het"
chr.data.lr.all = rbind(chr.data.lr.het, chr.data.lr.loh)
chr.data.lr.all$type <- ifelse(chr.data.lr.all$type == "Het", "No", "Yes")


https://urldefense.proofpoint.com/v2/url?u=http-3A__stat.lr&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=WV8Pnv3BPfAhAw6DReB65G63gatdIwCsUJyyNoCyU6Y&e=  <- as.data.frame(chr.data.lr.all) %>%
  wilcox_test(lr ~ type) %>%
  add_significance() %>% add_xy_position(x = "Response")
stat.lr$xmax = 2
stat.lr$xmin = 1
# stat.lr$y.position = 0
p <- chr.data.lr.all %>% ggplot(aes(x = type, y = lr))+
  geom_boxplot(width = 0.4)+
  geom_quasirandom(color = "#b90b59", varwidth = TRUE, size = 4, width = 0.21)+
  labs(y = "Chr17 Average Lr", x = "Chr17-LOH")+
  geom_hline(yintercept = 0.04613389, linetype = "dotted")+
  geom_hline(yintercept = -0.11, linetype = "dashed")+theme_pubr()+
  stat_pvalue_manual(https://urldefense.proofpoint.com/v2/url?u=http-3A__stat.lr&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=WV8Pnv3BPfAhAw6DReB65G63gatdIwCsUJyyNoCyU6Y&e= , label = "P = {p}",vjust = -0.01, bracket.nudge.y = -0.01, size = 7)+
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 25))+
  scale_y_continuous(breaks = seq(-0.5,0.45,0.1))
p
```


