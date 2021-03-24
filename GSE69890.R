###Ddownload the expression matrix and phenotype information files
phenotype=read.table('~/repo/phenotypes_euro.txt',sep = ',',header = T)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF/AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR/SR')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF/SR')
GSE69890=read.table('~/repo/GSE69890.txt',sep = '\t',header = T)


#annotate to keep ENsemble ID that have a official symbol name and belongs to the lncRNA and mRNA category as annotated from ENSEMBLE 99,
#and remove rows that have duplicated gene name
annotation=read.table("~/repo/annotation.txt",sep = '\t',header=T)
GSE69890$gene_id=rownames(GSE69890)
GSE69890=merge(GSE69890,annotation,by='gene_id')
GSE69890=GSE69890[GSE69890$gene_biotype%in%c('lncRNA','protein_coding'),]
row.names(GSE69890)=GSE69890$gene_id
GSE69890$gene_id=NULL
GSE69890=GSE69890[!duplicated(GSE69890$gene_name),]

#removing all features that have a count of less than  10 in more than 80% of the samples
GSE69890$num=apply(GSE69890[,1:235], 1, function(x) length(x[x>=10])>=235*0.2)
GSE69890=GSE69890[GSE69890$num=='TRUE',]

##DESeq2 for VST normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(GSE69890[,(1:235)]),
                              colData = phenotype,
                              design = ~AtrialRhythm)
rld <- varianceStabilizingTransformation(dds,blind = F)
vstMat <- assay(rld)
matrix=as.data.frame(vstMat)
matrix$gene_id=rownames(matrix)
matrix=merge(matrix,annotation,by='gene_id')
row.names(matrix)=matrix$gene_id
matrix$gene_id=NULL

###identifie SVs
library(sva)
mod = model.matrix(~-1+AtrialRhythm+sexFemale,data = phenotype)
mod0 = model.matrix(~sexFemale,data=phenotype)
n.sv = num.sv(matrix[,1:235],mod,method = 'leek')
svobj= sva(as.matrix(matrix[,1:235]),mod,mod0,n.sv=n.sv)

#regress out SVs and Sex
cleaningY = function(y, mod, svaobj) {
  X=cbind(mod,svaobj$sv)
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  P=ncol(mod)-1
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
dat.adjusted=cleaningY(as.matrix(matrix[,1:235]),mod,svobj)


#cleaned data
rt=data.frame(dat.adjusted)
rt$type=matrix$gene_biotype
rt$name=matrix$gene_name
write.table(rt,file = '~repo/rt.txt',sep = '\t')
