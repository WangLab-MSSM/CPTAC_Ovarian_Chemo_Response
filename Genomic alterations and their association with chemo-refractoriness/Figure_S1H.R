## Figure S1H:
```{r}
## this step needs raw data downloaded from dbgap and ran CNVEX using default setting
meta         <- fread("Meta data")
cnvex.purity <- fread("Purity ploidy measures file") ## we want purity values from this file

colnames(meta)[1] = "patient"
colnames(cnvex.purity) <- c("sample", "cand", "CNVEX.ploidy", "CNVEX.purity")

meta <- meta[patient %in% chr17.info$case]
meta <- merge.data.table(meta, cnvex.purity, all.x = T, by.x = "patient", by.y = "sample")

#### start the heatmap
oc.age = meta[, .(id=patient, age=PatientAge)]
oc.age = meta[, .(id=patient, NeoAdjuvant=NeoAdjuvant)]

digest.fns <- list.files("Reading all cnvex digest files",full.names = T, recursive = TRUE, pattern = "*-somatic-digest.rds")
fns <- data.table(digest.fn = digest.fns)
fns[, ids := str_replace(basename(digest.fn), "-somatic-digest.rds", "")]
fns[, patient := tstrsplit(ids,"-")[1]]
segs <- foreach(i = 1:nrow(fns)) %dopar% {
  dig <- readRDS(fns$digest.fn[i])
  seg = cnvex:::segOut(dig)
  cov <- dig$tile
  olap = findOverlaps((seg), (cov), select = "first")
  # seg$arm = cov$arm[olap]
  
  seg$sample <- fns$ids[i]
  seg$patient <- fns$patient[i]
  
  ## arm level
  return(as.data.table(seg))
}
segs <- rbindlist(segs)
.samplePd <- function(id,digest,tile.width) {
  ht.grid <- as.data.table(tileGenome(gobj$seqi, tilewidth=tile.width, cut.last.tile.in.chrom=TRUE)) ## heatmap grid
  colnames(ht.grid)[1] <- "chr"
  
  # segs = as.data.table(segOut(digest$seg, digest$fit))
  segs = as.data.table(segOut(digest))
  C.median = weighted.median(as.numeric(segs$C),as.numeric(segs$width),na.rm = TRUE)
  
  ## add events
  segs[,event.raw := fcase(
    C < 2, "Del",
    C > 2, "Gain",
    C == 2, "NoEvent"
  )]
  segs[,event.relative := fcase(
    C < C.median, "Del",
    C > C.median, "Gain",
    C == C.median, "NoEvent"
  )]
  segs[,loh := (K == 0),]
  segs[is.na(loh), loh := FALSE]
  
  ## add arm to segs
  cov <- as.data.table(digest$tile)
  colnames(cov)[1] <- "chr"
  colnames(segs)[1] <- "chr"
  segs <- .addCovToSeg(segs,cov)
  
  ## integrate with heatmap grid
  olap <- .dt.findOverlap(ht.grid,segs)
  ht.grid <- cbind(ht.grid, segs[olap,c("seg","C","K","hq","arm","event.raw","event.relative","loh")])
  ht.grid[, ':='(purity=digest$purity, ploidy=digest$ploidy,C.median = C.median, case=id)]
  
  ## TODO:info about arm level plot added here in a new object
  
  return(list(ht.pd = ht.grid))
}
.cohortPd <- function(add.list) {
  ## get heatmap plot data
  ht.pds = foreach(i = 1:length(add.list[["ids"]])) %dopar% {
    id <- add.list[["ids"]][i]
    digest <- readRDS(add.list[["digest.fns"]][i])
    ht.pd <- .samplePd(id,digest,tile.width)[["ht.pd"]]
    return(ht.pd)
  }
  ht.pd <- rbindlist(ht.pds)
  
  return(ht.pd)
}
tile.width = 1000000
## convert to drawable matrix
D.segs = segs[,.(D=weighted.mean(as.numeric(C),as.integer(nlr), na.rm=TRUE)),.(sample, patient)]
D.segs = merge.data.table(D.segs, meta, all.x = T, by = "patient")

setkey(D.segs, Response, patient)
D.segs = D.segs[order(Response, -D),]


## plot all
add.list <- .getAddresses(scratch,object.type = "digest")
add.list$ids <- str_replace(add.list$ids, "-wgs", "")
## make heatmap matrix
ht.pd <- .cohortPd(add.list) ## add tile width as a input here
ht.pd$chr = factor(str_replace(ht.pd$chr, "chr", ""), levels = str_replace(levels(ht.pd$chr), "chr",""))
ht.pd[, ids := str_replace(case, "-wgs", "")]
ht.pd[, patient := ids]
ht.pd$case = ht.pd$ids
ht.pd = ht.pd[!(chr %in% c("X", "Y")),]
https://urldefense.proofpoint.com/v2/url?u=http-3A__ht.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=G2t4xuBswFBcvJUUMtfsX8uvr9zJwZBHlgW_JVyAuQc&e=  <- .makeHeatmapMx(ht.pd,TRUE,TRUE) ## you can change this function case->ids to make name of each row better
# order
# ord <- match(D.segs$sample, rownames(https://urldefense.proofpoint.com/v2/url?u=http-3A__ht.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=G2t4xuBswFBcvJUUMtfsX8uvr9zJwZBHlgW_JVyAuQc&e= ))
ord <- sapply(1:nrow(D.segs), function(i) {
  idx = which(rownames(https://urldefense.proofpoint.com/v2/url?u=http-3A__ht.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=G2t4xuBswFBcvJUUMtfsX8uvr9zJwZBHlgW_JVyAuQc&e= ) == D.segs$patient[i])
  return(idx)
})

https://urldefense.proofpoint.com/v2/url?u=http-3A__ht.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=G2t4xuBswFBcvJUUMtfsX8uvr9zJwZBHlgW_JVyAuQc&e=  = ht.mx[ord,]

patient.split = do.call(rbind,strsplit(rownames(https://urldefense.proofpoint.com/v2/url?u=http-3A__ht.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=G2t4xuBswFBcvJUUMtfsX8uvr9zJwZBHlgW_JVyAuQc&e= ), "\\."))[,1]
# row.gaps = rep(0.06, length(patient.split))
row.gaps = c(rep(0.06, sum(D.segs$Response == "Refractory", na.rm=T)-1), 1.6, rep(0.06, sum(D.segs$Response == "Sensitive", na.rm=T)-1), 1.6, rep(0.06, sum(is.na(D.segs$Response))))

row.names <- rep("", length(row.gaps))
col_fun = colorRamp2(c(-9999,-3, -1, 0, 1, 5), c("#00ff4c","blue","#9382ff", "#e8e8e8", "#ff7675", "#c10002"))

ht.mx1 = https://urldefense.proofpoint.com/v2/url?u=http-3A__ht.mx&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-ALLUrcfR-4CCQkZVKC8w3o&r=APjqYeNFvbRrYlTnEjhIHxhnWWF3hg-ezYrOaGJn2wc&m=-6Fwzhu1uRqg_ihu_vLp9DrM0TfZZsn05vg5uQN90yR9wehnXrICPRldvVY8BV-I&s=G2t4xuBswFBcvJUUMtfsX8uvr9zJwZBHlgW_JVyAuQc&e= 
  patient.split1 = patient.split
row.gaps1 = row.gaps
D.segs[, Location := fcase(Location == "OV", "Ovary",
                           Location == "OM", "Omentum",
                           Location == "Mix", "Mix",
                           default = "Others")]
row_ha1 = rowAnnotation(Ploidy = D.segs$D, WGii = D.segs$wgii, nTAi = D.segs$tAi5, 
                        e1 = anno_empty(border = FALSE, width = unit(1, "mm")), ## instability measures
                        
                        Purity = D.segs$CNVEX.purity, Location = D.segs$Location,
                        e3 = anno_empty(border = FALSE, width = unit(1, "mm")), 
                        `Immune Factor` = D.segs$Immune_Factor, `Stroma Total Factor` = D.segs$StromaTotal_Factor ,
                        `Epithelial Factor` = D.segs$Epithelial_Factor, `Predicted Score`=D.segs$PredScore, `TCGA Subtype`= as.character(D.segs$TCGASubtype), 
                        e4 = anno_empty(border = FALSE, width = unit(1, "mm")), 
                        Response = D.segs$Response,
                        
                        
                        col = list(Ploidy   = colorRamp2(c(1.5,6), c("#e8e8e8", "#ff0800")),
                                   WGii     = colorRamp2(c(0,1), c("#e8e8e8", "#ff0800")),
                                   nTAi     = colorRamp2(c(10,44), c("#e8e8e8", "#ff0800")),
                                   
                                   Purity   = colorRamp2(c(0,1), c("#e8e8e8", "#26d945")),
                                   Location = c("Ovary" = "#26d945","Omentum" = "#0083ff", "Mix" = "#cc3395", "Others" = "#3b57c4"),
                                   
                                   `Immune Factor` = colorRamp2(c(0,1), c("#e8e8e8", "#66a100")),
                                   `Stroma Total Factor` = colorRamp2(c(0,1), c("#e8e8e8", "#66a100")),
                                   `Epithelial Factor` = colorRamp2(c(0,1), c("#e8e8e8", "#66a100")),
                                   `Predicted Score` = colorRamp2(c(0,1), c("#e8e8e8", "#66a100")),
                                   `TCGA Subtype` = c("DIF" = "#ff0028","IMR" = "#1803fc", "MES" = "#0083ff", "PRO" = "#622ed1"),
                                   
                                   Response = c("Sensitive"="#ff6d00", "Refractory"="#0083ff")
                                   
                        ),
                        show_legend=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold")
)
row.names1 = (row.names)

ht = ComplexHeatmap::Heatmap(ht.mx1, name = "CopyNumber",
                             cluster_columns = FALSE,cluster_rows = F, show_row_dend = F,show_heatmap_legend = F,
                             row_names_side = "left",
                             row_names_gp = gpar(fontsize = 7, fontface = "bold"),
                             row_title_gp = gpar(fontsize = 7, fontface = "bold"),
                             row_title = row.names1,
                             row_title_rot = 0,
                             column_names_gp = gpar(fontsize = 7, fontface = "bold"),
                             row_gap = unit(row.gaps1, "mm"),
                             show_row_names = FALSE,
                             row_split = factor(rownames(ht.mx1),levels = unique(rownames(ht.mx1))),
                             left_annotation = row_ha1,
                             column_split = ht.pd[case == ht.pd[1,case],chr],
                             column_title_gp = gpar(fontsize = 7, fontface = "bold"),
                             column_gap = unit(0.5, "mm"), 
                             heatmap_legend_param = list(),
                             col = col_fun,
                             na_col = "black",
                             use_raster = TRUE, raster_quality = 8)

```


