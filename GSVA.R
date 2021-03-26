exprs=read.table("~/repo/rt.txt",sep = "\t",header=T)
row.names(exprs)=exprs$name
exprs=exprs[,1:235]

library(GSVAdata)
library(Biobase)
library(RColorBrewer)
library(GSVA)
library(ggcorrplot)
library(ggthemes)
#download the gmt file
gmt_file="~/repo/c5.go.bp.v7.2.symbols.gmt"
geneset <- getGmt(gmt_file) 
GOBP=gsva(as.matrix(exprs),geneset,min.sz=10, max.sz=500, verbose=TRUE)
GOBP=as.data.frame(GOBP)

#phenotype
phenotype=read.table('~/repo/phenotypes_euro.txt',sep = ',',header = T)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF/AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR/SR')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF/SR')

#LINC00964 
gene=as.data.frame(t(exprs['LINC00964',]))
B1=c();B2=c();B3=c()
P1=c();P2=c();P3=c()
model=model.matrix(~AtrialRhythm,data = phenotype)
model=cbind(model,gene[,1])
for (i in 1:dim(GOBP)[1]) {
  phenotype$BP=as.numeric(GOBP[i,])
  fit=lm(BP~model,data = phenotype)
  fit=summary(fit)
  B1=c(B1,fit$coefficients[2,1])
  P1=c(P1,fit$coefficients[2,4])
  B2=c(B2,fit$coefficients[3,1])
  P2=c(P2,fit$coefficients[3,4])
  B3=c(B3,fit$coefficients[4,1])
  P3=c(P3,fit$coefficients[4,4])
}
BP=data.frame(B1,P1,B2,P2,B3,P3)
BP$BP=row.names(GOBP)
BP=BP[order(BP$B3,decreasing = T),]
write.table(BP,file='LINC00964_BP.txt',sep = '\t')

#MIAT 
gene=as.data.frame(t(exprs['MIAT',]))
B1=c();B2=c();B3=c()
P1=c();P2=c();P3=c()
model=model.matrix(~AtrialRhythm,data = phenotype)
model=cbind(model,gene[,1])
for (i in 1:dim(GOBP)[1]) {
  phenotype$BP=as.numeric(GOBP[i,])
  fit=lm(BP~model,data = phenotype)
  fit=summary(fit)
  B1=c(B1,fit$coefficients[2,1])
  P1=c(P1,fit$coefficients[2,4])
  B2=c(B2,fit$coefficients[3,1])
  P2=c(P2,fit$coefficients[3,4])
  B3=c(B3,fit$coefficients[4,1])
  P3=c(P3,fit$coefficients[4,4])
}
BP=data.frame(B1,P1,B2,P2,B3,P3)
BP$BP=row.names(GOBP)
BP=BP[order(BP$B3,decreasing = T),]
write.table(BP,file='MIAT_BP.txt',sep = '\t')











