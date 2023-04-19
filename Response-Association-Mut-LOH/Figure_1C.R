## Figure 1-C
```{r}
library(vcd)
library(gridGraphics)
library(gridExtra)
https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e=  <- fread("CHR17 LOH data") ## using this list of samples
mutations.relax <- fread("Mutation list for BRCA1 and BRCA2")

chr17.info$brca_mut <- chr17.info$case %in% mutations.relax[SYMBOL %in% c("BRCA1", "BRCA2")]$dna_id
chr17.info$brca_mut <- ifelse(chr17.info$brca_mut == TRUE, "Mut", "WT")

vcd::mosaic(Response ~ brca_mut, data = https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e= , gp=shading_max, split_vertical=T)  # brca

```



