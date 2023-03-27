## Figure 2-A pre-process:
```{r}
library(data.table)
library(doParallel)

rna.oc       <- fread("Normalized RNA expression data")
prot.oc      <- fread("Normalized RNA expression data")
prot.oc      <- fread("Chr17 LOH data")
meta         <- fread("Meta data")
https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e=      <- fread("Gene level Copy number values with median CN per sample moved to 0")
cnvex.purity <- fread("Purity ploidy measures file") ## we want purity values from this file



colnames(meta)[1] = "patient"
colnames(cnvex.purity) <- c("sample", "cand", "CNVEX.ploidy", "CNVEX.purity")

meta <- meta[patient %in% chr17.info$case]
meta <- merge.data.table(meta, cnvex.purity, all.x = T, by.x = "patient", by.y = "sample")

https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e=  <- na.omit(https://urldefense.proofpoint.com/v2/url?u=http-3A__cnv.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=Fextg7OOcxoUWGQ7WA5EypY-6U6AFHR29kx0SkIsqkM&e= )
## perform the analysis for RNA
https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e=  <- copy(rna.oc)

## getting rid of low-expressed genes
data.select.rowmean <- rowMeans(data.select[,-1])
thr                 <- quantile(data.select.rowmean, 0.3)
row.index           <- which(data.select.rowmean > thr)
https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e=  <- data.select[row.index, ]
doParallel::registerDoParallel(30)

## calculate association for all the genes
doParallel::registerDoParallel(30)
rna.cnv.association <- foreach(i = 1:nrow(https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e= )) %dopar% {
  
  ## expression
  row  <- (data.select[i, ])
  gene <- row$gene
  row  <- t(row[,-1])
  
  ## cnv
  cnv.gene <- cnv.data[gene_name == gene]
  if(nrow(cnv.gene) == 0) {
    return(NULL)
  }
  cnv.gene <- t(cnv.gene[, -(1:4)]) ## should be only numerical C values per sample
  cnv.gene <- data.table(dna = rownames(cnv.gene), cnv = cnv.gene[,1])
  cnv.gene <- merge.data.table(cnv.gene, meta[, c("patient", "Sample", "CNVEX.purity", "Cluster", "Response", "SampleAge")], by.x = "dna", by.y = "patient")
  
  ## combine
  https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- data.table(sample = rownames(row), exp = row[,1])
  https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , cnv.gene, by.x="sample", by.y = "Sample")
  
  
  ## apply purity as multiplicative
  reg.data$cnv.dose <- reg.data$cnv * reg.data$CNVEX.purity
  https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=           <- na.omit(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )
  reg.data$Response <- factor(reg.data$Response, levels = c("Refractory", "Sensitive"))
  
  ## correlation
  reg.data.split <- split(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , reg.data$Response)
  cor.r <- cor.test(reg.data.split[["Refractory"]]$exp, reg.data.split[["Refractory"]]$cnv.dose)
  cor.s <- cor.test(reg.data.split[["Sensitive"]]$exp, reg.data.split[["Sensitive"]]$cnv.dose)
  return(data.table(
    gene = gene, 
    r.cor.coef = cor.r$estimate , r.cor.pval = cor.r$p.value,
    s.cor.coef = cor.s$estimate , s.cor.pval = cor.s$p.value)
  )
}
rna.cnv.association <- rbindlist(rna.cnv.association)

rna.cnv.association$r.cor.qval <- p.adjust(rna.cnv.association$r.cor.pval, method = "BH")
rna.cnv.association$s.cor.qval <- p.adjust(rna.cnv.association$s.cor.pval, method = "BH")

###################################################################################################3
###################################################################################################3
###################################################################################################3
###########################  PROTEIN VS RNA
###################################################################################################3
###################################################################################################3
###################################################################################################3
https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e=   <- copy(prot.oc)
data.select2 <- copy(rna.oc)

## fix the boundry of genes
selected.genes     <- intersect(data.select$Index, data.select2$gene)
selected.samples   <- c(intersect(colnames(data.select2), colnames(https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e= )))
selected.samples.r <- c("gene",selected.samples)
selected.samples.p <- c("Index", selected.samples)
https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e=         <- data.select[Index %in% selected.genes]
data.select2       <- data.select2[gene %in% selected.genes]
https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e=         <- data.select[, ..selected.samples.p]
data.select2       <- data.select2[, ..selected.samples.r]

doParallel::registerDoParallel(30)
rna.prot.association <- foreach(i = 1:nrow(https://urldefense.proofpoint.com/v2/url?u=http-3A__data.select&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=ZMUshc3IiDbKL_3mnxwUEav6grgG-Zjhi1YwWeVPQGU&e= )) %dopar% {
  
  ## protein
  row  <- (data.select[i, ])
  https://urldefense.proofpoint.com/v2/url?u=http-3A__gene.name&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=h0h7D-24qgisqlJIvkdY8gsHKaDkNXskTmQ3pk8Q_Wo&e=  <- row$Index
  row  <- t(row[,-1])
  
  ## rna
  rna.gene <- data.select2[gene == https://urldefense.proofpoint.com/v2/url?u=http-3A__gene.name&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=h0h7D-24qgisqlJIvkdY8gsHKaDkNXskTmQ3pk8Q_Wo&e= ]
  if(nrow(rna.gene) == 0) {
    return(NULL)
  }
  rna.gene <- t(rna.gene[, -(1)])
  rna.gene <- data.table(Sample = rownames(rna.gene), rna = as.numeric(rna.gene[,1]))
  rna.gene <- merge.data.table(rna.gene, meta[, c("patient", "Sample", "CNVEX.purity", "Cluster", "Response", "SampleAge")], by.x = "Sample", by.y = "Sample")
  
  ## combine
  https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- data.table(sample = rownames(row), prot = as.numeric(row[,1]))
  https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=  <- merge.data.table(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , rna.gene, by.x="sample", by.y = "Sample")
  
  https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e=           <- na.omit(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= )
  reg.data$Response <- factor(reg.data$Response, levels = c("Refractory", "Sensitive"))
  
  reg.data.split <- split(https://urldefense.proofpoint.com/v2/url?u=http-3A__reg.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=qiokNcTfccdQbfL7ECGoeAPMT7JWWL3fd3SErmUt0xQ&e= , reg.data$Response)
  
  ## correlation:
  cor.r <- cor.test(reg.data.split[["Refractory"]]$prot, reg.data.split[["Refractory"]]$rna)
  cor.s <- cor.test(reg.data.split[["Sensitive"]]$prot, reg.data.split[["Sensitive"]]$rna)
  
  return(data.table(
    gene = https://urldefense.proofpoint.com/v2/url?u=http-3A__gene.name&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=h0h7D-24qgisqlJIvkdY8gsHKaDkNXskTmQ3pk8Q_Wo&e= , 
    r.cor.coef = cor.r$estimate , r.cor.pval = cor.r$p.value,
    s.cor.coef = cor.s$estimate , s.cor.pval = cor.s$p.value
  ))
}
rna.prot.association <- rbindlist(rna.prot.association)

rna.prot.association$r.cor.qval <- p.adjust(rna.prot.association$r.cor.pval, method = "BH")
rna.prot.association$s.cor.qval <- p.adjust(rna.prot.association$s.cor.pval, method = "BH")
```

## Figure 2-A
```{r}
library(tidyverse)
library(rstatix)
library(ggpubr)

boxplot.data.rnaPRot <- data.table(type = c(rep("Ref", nrow(rna.prot.association)), rep("Sen", nrow(rna.prot.association))),
                                   cor  = c(rna.prot.association$r.cor.coef, rna.prot.association$s.cor.coef),
                                   qval = c(rna.prot.association$r.cor.qval, rna.prot.association$s.cor.qval),
                                   gene = c(rna.prot.association$gene, rna.prot.association$gene))

boxplot.data.cnvRna <- data.table(type = c(rep("Ref", nrow(rna.cnv.association)), rep("Sen", nrow(rna.cnv.association))),
                                  
                                  cor  = c(rna.cnv.association$r.cor.coef, rna.cnv.association$s.cor.coef),
                                  qval = c(rna.cnv.association$r.cor.qval, rna.cnv.association$s.cor.qval),
                                  gene = c(rna.cnv.association$gene, rna.cnv.association$gene))


boxplot.data.cnvRna <- boxplot.data.cnvRna[gene %in% boxplot.data.rnaPRot$gene]
boxplot.data.rnaPRot <- boxplot.data.rnaPRot[gene %in% boxplot.data.cnvRna$gene]

boxplot.data.cnvRna$class     <- "CNV-RNA"
boxplot.data.rnaPRot$class <- "RNA-Protein"

boxplot.data.cnvRna     <- boxplot.data.cnvRna[order(gene, type)]
boxplot.data.rnaPRot <-boxplot.data.rnaPRot[order(gene, type)]

boxplot.data.cnvRna <- boxplot.data.cnvRna[, c("type", "cor", "qval", "gene", "class")]
https://urldefense.proofpoint.com/v2/url?u=http-3A__boxplot.data&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=EwqM7818uTP74XgvgkqsdOcKx3D-ZK8zAiekX8JvTUg&e=  <- rbind(boxplot.data.cnvRna, boxplot.data.rnaPRot)


pd1 <- copy(boxplot.data.cnvRna)
pd1$prot.cor <- boxplot.data.rnaPRot$cor
pd1$pd <- pd1$cor - pd1$prot.cor

stat.cnvrna <- as.data.frame(boxplot.data[class == "CNV-RNA"]) %>%
  wilcox_test(cor ~ type) %>%
  add_significance() %>% add_xy_position(x = "type")
p1 = as.data.frame(boxplot.data[class == "CNV-RNA"]) %>% ggplot(aes(x = type, y = cor, color = type))+
  geom_boxplot()+scale_color_manual(values = c("#0073C2", "#E7B800"))+
  stat_pvalue_manual(stat.cnvrna, label = "P = {p}", vjust = 0, bracket.nudge.y = 0.05, size = 7)+
  theme_pubr()+
  labs(x = "Response", y = "CNV-RNA\nCorrelation")+
  theme(axis.title = element_text(size = 20), axis.text =  element_text(size = 18))+
  coord_cartesian(ylim = c(-0.5,1))

stat.rnaprotein <- as.data.frame(boxplot.data[class == "RNA-Protein"]) %>%
  wilcox_test(cor ~ type) %>%
  add_significance() %>% add_xy_position(x = "type")

p2 = as.data.frame(boxplot.data[class == "RNA-Protein"]) %>% ggplot(aes(x = type, y = cor, color = type))+
  geom_boxplot()+scale_color_manual(values = c("#0073C2", "#E7B800"))+
  stat_pvalue_manual(stat.rnaprotein, label = "P = {p}", vjust = 0, bracket.nudge.y = 0.05, size = 7)+
  theme_pubr()+
  labs(x = "Response", y = "CNV-Protein\nCorrelation")+
  theme(axis.title = element_text(size = 20), axis.text =  element_text(size = 18), 
        axis.text.y = element_blank())+
  coord_cartesian(ylim = c(-0.5,1))

p1 | p2

```



