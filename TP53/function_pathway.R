# 


# function for pathway enrichment analysis with p-values  -----------------------------------------------------------------------


logp_pathway1 = function(GeneID=reg.summary$Gene_ID, 
         EST=reg.summary$estimate, 
         Pvalue=reg.summary$p.value, 
         pathway="Hallmark") {
  if (pathway=="Hallmark") {
    path.database=GSA::GSA.read.gmt("Hallmark.gmt")
  }
  if (pathway== "DDR") {
    path.database=GSA::GSA.read.gmt("DDR.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway== "KEGG") {
    path.database=GSA::GSA.read.gmt("KEGG.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway== "DQ") {
    path.database= GSA::GSA.read.gmt("PTRC_DQ_list.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway== "Reactome") {
    path.database= GSA::GSA.read.gmt("Reactome.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway=="TGFb") {
    path.database= GSA::GSA.read.gmt("TGFb.gmt")
    path.database$geneset.names=paste0("TGFb_", path.database$geneset.names)
  }

  score= sign(EST) * (-log10(Pvalue))
  score.rank=rank(score)  
  path.p = path.direction = avglogFC = n.path=rep(NA, length(path.database$genesets))

  for (i in 1:length(path.database$genesets)){
    idx=which(GeneID %in% path.database$genesets[[i]])
    avglogFC[i]=mean(EST[idx])
    n.path[i]=length(idx)

    score.in <- score.rank[idx]
    score.out <- score.rank[-idx]
    if (length(score.in>=2)) {
      test=wilcox.test(score.in, score.out, alternative = "two.sided")
      path.p[i]=test$p.value
      path.direction[i]=ifelse (mean(score.in) > mean(score.out), 1, -1)
    } else {path.p[i]=path.direction[i]=NA}
   
    print(i)
  }
  fdr.BY=p.adjust(path.p, method="BY")
  # signed.log.fdr=path.direction* (-log10(fdr.BY))
  return(data.frame(geneset.names=path.database$geneset.names, 
               path.direction=path.direction,
               path.p=path.p, 
               fdr.BY=fdr.BY,
               avglogFC=avglogFC,
               n.path=n.path))
} 

# function for pathway enrichment analysis with p-values  -----------------------------------------------------------------------

logp_pathway_multiple = function(GeneID=reg.summary$Gene_ID, 
                         EST=reg.summary$estimate, 
                         Pvalue=reg.summary$p.value,
                         pathways=c("Hallmark", "KEGG", "Reactome", "DDR", "DQ", "TGFb")) {
  lp=length(pathways)
  p.all=NULL
  for (i in 1:lp) {
    p1=logp_pathway1(GeneID=GeneID, 
                    EST=EST, 
                    Pvalue=Pvalue, 
                    pathway=pathways[i])
    p.all=rbind(p.all, p1)
  }
  return(p.all)
}

# function for pathway enrichment analysis with p-values  -----------------------------------------------------------------------


binary_pathway1 = function(YesList=Genes_R_cascade, NoList=Genes_R_not, 
                         pathway="Hallmark") {
  if (pathway=="Hallmark") {
    path.database=GSA::GSA.read.gmt("Hallmark.gmt")
  }
  if (pathway== "DDR") {
    path.database=GSA::GSA.read.gmt("DDR.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway== "KEGG") {
    path.database=GSA::GSA.read.gmt("KEGG.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway== "DQ") {
    path.database= GSA::GSA.read.gmt("PTRC_DQ_list.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway== "Reactome") {
    path.database= GSA::GSA.read.gmt("Reactome.gmt")
    path.database$geneset.names=paste0(path.database$geneset.descriptions, "_", path.database$geneset.names)
  }
  if (pathway=="TGFb") {
    path.database= GSA::GSA.read.gmt("TGFb.gmt")
    path.database$geneset.names=paste0("TGFb_", path.database$geneset.names)
  }

  n_path = n_path_data= prop.yes.in.path=or= 
    prop.yes.outside.path=fisher.pvalue =  rep(NA, length(path.database$genesets))

  for (i in 1:length(path.database$genesets)){
    path_gene_once <- setdiff(path.database$genesets[[i]], "")
    n_path[i] <- length(path_gene_once) 
    n_path_data[i] <-length(intersect(path_gene_once, c(YesList, NoList)) )

    table <- matrix(c(sum(YesList %in% path_gene_once),
                              sum(!YesList %in%  path_gene_once),
                              sum(NoList %in% path_gene_once),
                              sum(!NoList %in% path_gene_once)), ncol=2)
    prop.yes.in.path[i] = table[1,1]/(sum(table[1,]))
    prop.yes.outside.path[i] = table[2,1]/(sum(table[2,]))    
    or.table=oddsratio.wald(table)
    or[i]=or.table$measure[2,1]
    fisher.pvalue[i]=or.table$p.value[2,2]
    
  }
  fdr.BY=p.adjust(fisher.pvalue, method="BY")
  result_combined <- tibble(
                              geneset.names = path.database$geneset.names,
                              n_path=n_path,
                              n_path_data=n_path_data,
                              prop.yes.in.path=prop.yes.in.path,
                            prop.yes.outside.path=prop.yes.outside.path,
                              or = or,
                              fisher.pvalue=fisher.pvalue, fdr.BY=fdr.BY)
  
  return(result_combined)
  
}


# function for pathway enrichment analysis with p-values  -----------------------------------------------------------------------

binary_pathway_multiple = function(YesList=Genes_R_cascade, NoList=Genes_R_not,
                                   pathways=c("Hallmark", "KEGG", "Reactome", "DDR", "DQ", "TGFb"))  {
  lp=length(pathways)
  p.all=NULL
  for (i in 1:lp) {
    p1=binary_pathway1(YesList=YesList, 
                       NoList=NoList, 
                    pathway=pathways[i])
    p.all=rbind(p.all, p1)
  }
  return(p.all)

}

