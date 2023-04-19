## Figure 1B
```{r}
https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e=  <- fread("chr17 loh data") ## using this list of samples
mutations.relax  <- fread("mutation file")
mutations.relax  <- mutations.relax[dna_id %in% chr17.info$case, c("dna_id","SYMBOL", "Consequence", "IMPACT")]

# mutations.relax[, c("Consequence", "IMPACT", "clinvar")]

## first the single consequences
mutations.relax[, Consequence := tstrsplit(Consequence, "&")[[1]]]
mutations.relax[, mut_type := ifelse(str_detect(Consequence, "frameshift"), "Frameshift", NA)]
mutations.relax[, mut_type := ifelse(str_detect(Consequence, "missense"), "Missense Mutation", mut_type)]
mutations.relax[, mut_type := ifelse(str_detect(Consequence, "splice") & IMPACT == "HIGH", "Splice site", mut_type)]
mutations.relax[, mut_type := ifelse(str_detect(Consequence, "stop_gained"), "Stop gained", mut_type)]
mutations.relax[, mut_type := ifelse(str_detect(Consequence, "inframe"), "Inframe InDels", mut_type)]
mutations.relax[, mut_type := ifelse(is.na(mut_type), "Other", mut_type)]

mutations.relax$IMPACT <- factor(mutations.relax$IMPACT, levels = c("MODERATE", "HIGH"), ordered = TRUE)
mutations.relax  <- mutations.relax[,.SD[which.max(IMPACT)], by = .(SYMBOL, dna_id)]
setkey(mutations.relax, "dna_id")
mutations.relax <- mutations.relax[J(dna_id = (chr17.info$case))]


mutations.relax.casted <- dcast.data.table(mutations.relax, SYMBOL ~ dna_id, value.var = "mut_type")
rnames <- mutations.relax.casted$SYMBOL
mutations.relax.casted <- as.matrix(mutations.relax.casted[,-1])
rownames(mutations.relax.casted) <- rnames
col <- c(Frameshift = "#E64B35FF",
         `Missense Mutation` = "#4DBBD5FF",
         `Splice site` = "#00A087FF",
         `Stop gained` = "#3C5488FF",
         `Inframe InDels` = "#F39B7FFF",
         Other = "#91D1C2FF")

# mutations.relax.casted <- mutations.relax.casted[!is.na(rownames(mutations.relax.casted)),]
mutations.relax.casted.sen <- mutations.relax.casted[,colnames(mutations.relax.casted) %in% chr17.info[Response == "Sensitive"]$case]
mutations.relax.casted.ref <- mutations.relax.casted[,colnames(mutations.relax.casted) %in% chr17.info[Response == "Refractory"]$case]


### sensitive
sen.tbl      <- as.data.table(t(mutations.relax.casted.sen))
count.mx.sen <- matrix(!is.na(sen.tbl), nrow = nrow(sen.tbl))
gene.freq.sen    <- data.table(Gene = colnames(sen.tbl), Freq = colMeans(count.mx.sen))
right.annot <- rowAnnotation(`Freq S` = gene.freq.sen$Freq,
                             col = list(`Freq S` = colorRamp2(c(0,0.2), c("white", "#26d945"))), show_legend = TRUE)


top.anno.tbl <- data.table(case = colnames(mutations.relax.casted.sen))
top.anno.tbl <- merge.data.table(top.anno.tbl, https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e= )

top_annotation_sen = HeatmapAnnotation(cbar = anno_oncoprint_barplot(axis_param = list(gp = gpar(fontsize = 12))),
                                       Chr17LOH = ifelse(top.anno.tbl$group == "LOH", "Yes", "No"),
                                       nTAI = top.anno.tbl$nTAI,
                                       Response = top.anno.tbl$Response,
                                       
                                       col = list(
                                         Chr17LOH = c("Yes" = "#0083ff", "No" = "#e8e8e8"),
                                         nTAI     = colorRamp2(c(0,44), c("#e8e8e8", "#ff0800")),
                                         Response = c("Sensitive"="#E7B800", "Refractory"="#0073C2")
                                       ),
                                       
                                       annotation_name_side = "left",
                                       
                                       annotation_name_gp = gpar(fontsize = 12),
                                       show_legend = FALSE
                                       
)

onco.sen <- oncoPrint(mutations.relax.casted.sen,
                      
                      alter_fun = list(
                        background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                                                    gp = gpar(fill = "#e8e8e8", col = "white")),
                        Frameshift = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                    gp = gpar(fill = col["Frameshift"], col = NA)),
                        `Missense Mutation` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                             gp = gpar(fill = col["Missense Mutation"], col = NA)),
                        `Splice site` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                       gp = gpar(fill = col["Splice site"], col = NA)),
                        `Stop gained` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                       gp = gpar(fill = col["Stop gained"], col = NA)),
                        `Inframe InDels` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                          gp = gpar(fill = col["Inframe InDels"], col = NA)),
                        Other = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                               gp = gpar(fill = col["Other"], col = NA))
                      ), col = col,
                      # pct_side = "right",
                      # row_names_side = "left",
                      remove_empty_columns = FALSE,
                      left_annotation =  rowAnnotation(
                        rbar = anno_oncoprint_barplot(
                          axis_param = list(direction = "reverse", gp = gpar(fontsize = 12)),
                        )),
                      right_annotation = right.annot,
                      top_annotation = top_annotation_sen,
                      row_names_gp =  gpar(fontsize = 12),
                      pct_gp =  gpar(fontsize = 12),
                      show_heatmap_legend = FALSE
)


## Refractory
ref.tbl      <- as.data.table(t(mutations.relax.casted.ref))
count.mx.ref <- matrix(!is.na(ref.tbl), nrow = nrow(ref.tbl))
gene.freq.ref    <- data.table(Gene = colnames(ref.tbl), Freq = colMeans(count.mx.ref))
left.annot <- rowAnnotation(`Freq R` = gene.freq.ref$Freq,
                            col = list(`Freq R` = colorRamp2(c(0,0.2), c("white", "#26d945"))), show_legend = TRUE)

top.anno.tbl.ref <- data.table(case = colnames(mutations.relax.casted.ref))
top.anno.tbl.ref <- merge.data.table(top.anno.tbl.ref, https://urldefense.proofpoint.com/v2/url?u=http-3A__chr17.info&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=OHHyKRa-9dVkqT1VPplBoTWToY2n7mY8PelYLgJuSjc&e= )

top_annotation_ref = HeatmapAnnotation(cbar = anno_oncoprint_barplot(axis_param = list(gp = gpar(fontsize = 12))),
                                       Chr17LOH = ifelse(top.anno.tbl.ref$group == "LOH", "Yes", "No"),
                                       nTAI = top.anno.tbl.ref$nTAI,
                                       Response = top.anno.tbl.ref$Response,
                                       
                                       col = list(
                                         Chr17LOH = c("Yes" = "#0083ff", "No" = "#e8e8e8"),
                                         nTAI     = colorRamp2(c(0,44), c("#e8e8e8", "#ff0800")),
                                         Response = c("Sensitive"="#E7B800", "Refractory"="#0073C2")
                                       ),
                                       
                                       annotation_name_side = "right",
                                       annotation_name_gp = gpar(fontsize = 12),
                                       show_legend = FALSE
)

onco.ref <- oncoPrint(mutations.relax.casted.ref,
                      alter_fun = list(
                        background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                                                    gp = gpar(fill = "#e8e8e8", col = "white")),
                        Frameshift = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                    gp = gpar(fill = col["Frameshift"], col = NA)),
                        `Missense Mutation` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                             gp = gpar(fill = col["Missense Mutation"], col = NA)),
                        `Splice site` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                       gp = gpar(fill = col["Splice site"], col = NA)),
                        `Stop gained` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                       gp = gpar(fill = col["Stop gained"], col = NA)),
                        `Inframe InDels` = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                                          gp = gpar(fill = col["Inframe InDels"], col = NA)),
                        Other = function(x, y, w, h) grid.rect(x, y, w*0.92, h*0.92, 
                                                               gp = gpar(fill = col["Other"], col = NA))
                      ), col = col,
                      pct_side = "right",
                      row_names_side = "left",
                      remove_empty_columns = FALSE,
                      left_annotation = left.annot,
                      top_annotation = top_annotation_ref,
                      right_annotation =  rowAnnotation(
                        rbar = anno_oncoprint_barplot(
                          axis_param = list(side = "bottom", gp = gpar(fontsize = 12))
                        )
                      ),
                      row_names_gp =  gpar(fontsize = 12),
                      pct_gp =  gpar(fontsize = 12),
                      show_heatmap_legend = FALSE
)


ComplexHeatmap::draw(onco.sen + onco.ref )
```
