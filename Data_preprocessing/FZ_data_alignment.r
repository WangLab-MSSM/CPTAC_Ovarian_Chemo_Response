
FD.bc = read.delim2('FD_GLBLprot_MI_FDbridge_Abund_20201002_Batch_v02.tsv',sep ='\t',header = T,row.names = 1,check.names = F)
FD.im = read.delim2('FD_GLBLprot_MI_FDbridge_Abund_20201002_Imput_v02.tsv',sep ='\t',header = T,row.names = 1,check.names = F)

FZ.bc = read.delim2('FZ_GLBLprot_MI_FZbridge_Abund_20201002_Batch_v02.tsv',sep ='\t',header = T,row.names = 1,check.names = F)
FZ.im = read.delim2('FZ_GLBLprot_MI_FZbridge_Abund_20201002_Imput_v02.tsv',sep ='\t',header = T,row.names = 1,check.names = F)

fd.rep = as.data.frame(FD.bc[,-(1:8)][,match(intersect(substr(colnames(FD.bc)[-(1:8)],8,12),substr(colnames(FZ.bc)[-(1:8)],8,12)),substr(colnames(FD.bc)[-(1:8)],8,12))])
fz.rep = as.data.frame(FZ.bc[,-(1:8)][,match(intersect(substr(colnames(FD.bc)[-(1:8)],8,12),substr(colnames(FZ.bc)[-(1:8)],8,12)),substr(colnames(FZ.bc)[-(1:8)],8,12))])

for(i in 1:dim(fd.rep)[2])
{fd.rep[,i] = as.numeric(fd.rep[,i])}

for(i in 1:dim(fz.rep)[2])
{fz.rep[,i] = as.numeric(fz.rep[,i])}

p.list = intersect(rownames(fd.rep),rownames(fz.rep))

m0.d = apply(fd.rep[p.list,],1,mean,na.rm = T)
s0.d = apply(fd.rep[p.list,],1,sd,na.rm = T)
m0.z = apply(fz.rep[p.list,],1,mean,na.rm = T)
s0.z = apply(fz.rep[p.list,],1,sd,na.rm = T)

coef.0 = cbind(m0.d-m0.z*s0.d/s0.z,s0.d/s0.z)

rownames(coef.0) = p.list

FZ.0 = matrix(as.numeric(as.matrix(FZ.im[rownames(coef.0),-(1:8)])),nrow(coef.0))*coef.0[2,]+coef.0[1,]

colnames(FZ.0) = colnames(FZ.im)[-(1:8)]
rownames(FZ.0) = p.list

FZ.out = data.frame(FZ.im[rownames(coef.0),(1:8)],FZ.0)

rownames(FZ.out) = rownames(coef.0)

write.table(FZ.out,'FZ_GLBLprot_MI_FZbridge_Abund_20201002_Imput_v02_AlignedToFD.tsv',sep = '\t')


