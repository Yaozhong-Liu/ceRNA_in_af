###GSE69890
phenotype=read.table('G:/生信分析/ceRNA/2020.12/GSE69890_RAW/phenotypes_euro.txt',sep = ',',header = T,row.names = 1)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF/AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR/SR')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF/SR')
table(phenotype$AtrialRhythm)
#library(tidyverse)
#phenotype=separate(phenotype,col=AtrialRhythm,into=c('History_of_AF','AF_at_surgery'),sep = '/')
#a=read.table('clipboard',sep = ',',header = T)
#GSE69890=a[,colnames(a)%in%colnames(phenotype)]
#write.table(GSE69890,file = 'G:/生信分析/ceRNA/2020.12/GSE69890_RAW/GSE69890.txt',sep = '\t')
GSE69890=read.table('G:/生信分析/ceRNA/2020.12/GSE69890_RAW/GSE69890.txt',sep = '\t',header = T,row.names = 1)


###annotate to remove ENsemble ID that have a official symbol name and belongs to the lncRNA and mRNA category as annotated from ENSEMBLE 99,
###and remove rows that have duplicated gene name
annotation=read.table("G:/生信分析/ceRNA/Probe_annotation/annotation.txt",sep = '\t',header=T)
GSE69890$gene_id=rownames(GSE69890)
GSE69890=merge(GSE69890,annotation,by='gene_id')
GSE69890=GSE69890[GSE69890$gene_biotype%in%c('lncRNA','protein_coding'),]
row.names(GSE69890)=GSE69890$gene_id
GSE69890$gene_id=NULL
GSE69890=GSE69890[!duplicated(GSE69890$gene_name),]

#removing all features that have a count of less than  10 in more than 80% of the samples
GSE69890$num=apply(GSE69890[,1:235], 1, function(x) length(x[x>=10])>=235*0.2)
###Remove rows that have low counts <50
#GSE69890=GSE69890[GSE69890$sum>=50,]
GSE69890=GSE69890[GSE69890$num=='TRUE',]
table(GSE69890$gene_biotype)#5731 lncRNA 17157 protein_coding


###DESEQ2 FOR VST normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(GSE69890[,(1:235)]),
                              colData = phenotype,
                              design = ~AtrialRhythm)

rld <- varianceStabilizingTransformation(dds,blind = F)
vstMat <- assay(rld)
matrix=as.data.frame(vstMat)
#boxplot(matrix,outline=F)
write.csv(matrix,file = 'G:/生信分析/ceRNA/2020.12/GSE69890_RAW/GSE69890_matrix_all.csv')
######################################################################################################################################
#BiocManager::install('GenomeInfoDbData')
matrix=read.csv('G:/生信分析/ceRNA/2020.12/GSE69890_RAW/GSE69890_matrix_all.csv',header = T,row.names = 1)
phenotype=read.table('G:/生信分析/ceRNA/2020.12/GSE69890_RAW/phenotypes_euro.txt',sep = ',',header = T,row.names = 1)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF/AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR/SR')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF/SR')


annotation=read.table("G:/生信分析/ceRNA/Probe_annotation/annotation.txt",sep = '\t',header=T)
matrix$gene_id=rownames(matrix)
matrix=merge(matrix,annotation,by='gene_id')
row.names(matrix)=matrix$gene_id
matrix$gene_id=NULL
table(matrix$gene_biotype)

###identifie SVs
library(sva)
mod = model.matrix(~-1+AtrialRhythm+sexFemale,data = phenotype)
mod0 = model.matrix(~sexFemale,data=phenotype)
n.sv = num.sv(matrix[,1:235],mod,method = 'leek')
svobj= sva(as.matrix(matrix[,1:235]),mod,mod0,n.sv=n.sv)
plot(svobj$sv)


###for limma
design <- cbind(mod,svobj$sv)
colnames(design)=c('AFAF','AFSR','SRSR','SEX','SV1','SV2')
library(limma)
edata=matrix
edata=edata[,1:235]
fit = lmFit(edata,design)
contrast.matrix=makeContrasts(contrasts = c("AFAF-AFSR",'AFSR-SRSR'),levels = design)
fit1=contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit1)
allDiff1=topTable(fit2,adjust='fdr',coef="AFAF-AFSR",number=200000)
allDiff2=topTable(fit2,adjust='fdr',coef="AFSR-SRSR",number=200000)
write.table(allDiff1,file = 'G:/生信分析/ceRNA/2020.12/alldiff_persistence.txt',sep = '\t')
write.table(allDiff2,file = 'G:/生信分析/ceRNA/2020.12/alldiff_susceptible.txt',sep = '\t')




cleaningY = function(y, mod, svaobj) {
  X=cbind(mod,svaobj$sv)
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  P=ncol(mod)-1
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
###regress out SVs and Sex

dat.adjusted=cleaningY(as.matrix(matrix[,1:235]),mod,svobj)



###clean data with regression
#model=cbind(mod,svobj$sv)
#dat=as.matrix(matrix[,1:235])
#dat.adjusted = dat
#for (i in 1:dim(dat)[1]) {
#  dat.adjusted[i,] <- as.matrix(model)[,1:3]%*%c(lm(dat[i,]~model)$coefficients['(Intercept)'],
#                                                 lm(dat[i,]~model)$coefficients[3],
##                                                 lm(dat[i,]~model)$coefficients[4],
#                                                 lm(dat[i,]~model)$coefficients[5]
#                                               )[1:3]+lm(dat[i,] ~model)$residuals}

 
###Combat to adjust for Sex
#batch = phenotype$sexFemale
#modcombat = model.matrix(~AtrialRhythm, data=phenotype)
#dat.adjusted= ComBat(dat=as.matrix(matrix[,1:235]), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


##---------------WGCNA with cleaned data
rt=data.frame(dat.adjusted)
rt$type=matrix$gene_biotype
rt$name=matrix$gene_name
table(rt$type)
write.table(rt,file = 'G:/生信分析/ceRNA/2020.12/rt.txt',sep = '\t')

rt=read.table('G:/生信分析/ceRNA/2020.12/rt.txt',sep = '\t',header = T,row.names = 1)
#boxplot(rt[,1:235],outline=F)
rt$var=apply(rt[,1:235],1,var);rt=rt[order(rt$var,decreasing = T),]
#mrna=rt[rt$type=='protein_coding',]
#nonmrna=rt[!rt$type=='protein_coding',]
#datE=rbind(mrna[1:4000,],nonmrna[1:1000,])
#datE=rbind(mrna[1:round(length(row.names(mrna))/4),],nonmrna[1:round(length(row.names(nonmrna))/4),])

#datE=rt[1:round(nrow(rt)*0.3),]
datE=rt[1:5000,]
table(datE$type)




phenotype=read.table('G:/生信分析/ceRNA/2020.12/GSE69890_RAW/phenotypes_euro.txt',sep = ',',header = T,row.names = 1)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF/AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR/SR')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF/SR')

