annotation=read.table("~/repo/annotation.txt",sep = '\t',header=T)
colnames(annotation)=c('geneId','genesymbol','genebiotype')

#Tarbase8.0
library(dplyr)
Tarbase=read.table('~/repo/TarBase_v8_download.txt',sep = '\t',header = T)
Tarbase=Tarbase[Tarbase$species=='Homo sapiens',]
Tarbase=Tarbase[,1:3]
Tarbase=distinct(Tarbase)#422893
Tarbase=merge(Tarbase,annotation,by='geneId')#419328
lnc_Tar=Tarbase[Tarbase$genebiotype=='lncRNA',]#645 pairs
pro_Tar=Tarbase[Tarbase$genebiotype=='protein_coding',]#417860 pairs


#lncbase2.0 experimental validated
lncbase=read.csv('~/repo/LncBasev2_download.csv',header = F)
library(tidyverse)
lncbase=separate(data = lncbase, col = V1, into = c('geneId','geneName','mirna','species','cell_line','tissue','category','tmethod','positive_negative','direct_indirect','condition'), sep = "\t")
lncbase=lncbase[-1,]
lncbase=lncbase[lncbase$species=='Homo sapiens',]
lncbase=lncbase[,c(1:3)]
lncbase=distinct(lncbase)#62273
lncbase=merge(lncbase,annotation,by='geneId')#47404
lnc_lnc=lncbase[lncbase$genebiotype=='lncRNA',]#25633 pairs
pro_lnc=lncbase[lncbase$genebiotype=='protein_coding',]#906 pairs

#rewrite
Tarbase=rbind(pro_Tar,pro_lnc)
Tarbase=Tarbase[,1:3]
Tarbase=distinct(Tarbase)#418758 pairs

lncbase=rbind(lnc_lnc,lnc_Tar)
lncbase=lncbase[,1:3]
lncbase=distinct(lncbase)#26178 pairs




#ready fpr input
key=read.table("~/repo/key.txt",sep = '\t',header = T)
MEs=read.table("~/repo/MEs.txt",sep = '\t',header = T)
library(stringr)
color=unique(colnames(MEs))
color=str_remove(color,'ME')
GO=c()
REA=c()

library('tcltk')
pb <- tkProgressBar("进度","已完成 %", 0, 100)
#u=1:(length(color)-1)
u=c(1,5,16,18)#4 AF related modules
for (LYZ in u){
lyz=color[LYZ]
key1=key[key$moduleColors==lyz,]
table(key1$type)
mRNA=key1[key1$type=='protein_coding',]
lncRNA=key1[!key1$type=='protein_coding',]

#ready for lnc and mrna
lnc=lncbase[lncbase$geneId%in%row.names(lncRNA),]
mrna=Tarbase[Tarbase$geneId%in%row.names(mRNA),]

#lncRNA-miRNA-mRNA pairs
cerna=merge(lnc,mrna,by='mirna')
cerna=cerna[,c(2,1,4)]
colnames(cerna)=c('lncID','mirna','mrnaID')
cerna$mirna=as.character(cerna$mirna)
cerna$lncID=as.character(cerna$lncID)
cerna$mrnaID=as.character(cerna$mrnaID)


N=length(unique(unique(lncbase$mirna),unique(Tarbase$mirna)))#923
#hypergeometric test
library(stats)
pvalue=c()
for (i in 1:length(row.names(cerna))) {
  LNCRNA=cerna$lncID[[i]]
  MRNA=cerna$mrnaID[[i]]
  t=cerna[cerna$lncID==LNCRNA&cerna$mrnaID==MRNA,]
  I=lncbase[lncbase$geneId==LNCRNA,];I=unique(as.character(I$mirna))
  m=Tarbase[Tarbase$geneId==MRNA,];m=unique(as.character(m$mirna))
  p=phyper(length(unique(t$mirna))-1,length(m),N-length(m),length(I),lower.tail = F)
  pvalue=c(pvalue,p)
}
cerna$p=pvalue

#BH and write
cerna1=cerna[,c(1,3,4)]
cerna1=distinct(cerna1)
library(qvalue)
cerna1$adjp=p.adjust(cerna1$p, method = 'BH', n = length(cerna1$p))
cerna1$merge=paste(cerna1$lncID,cerna1$mrnaID)
cerna$merge=paste(cerna$lncID,cerna$mrnaID)
final=merge(cerna,cerna1,by='merge')
final=final[,c(2,3,4,8,9)]
colnames(final)=c('lncID','mirna','mrnaID','p','adjusted p')
write.table(final,file = paste0('~/repo/all_triple_pair_',lyz,'.txt'),sep = '\t',row.names = F)


#lncRNA-mRNA pair
pair=final[,c(1,3,4,5,6)]
colnames(pair)=c('node1','node2','p','adjustedp','type')
allpair=pair
allpair=distinct(allpair)
symbolnode1=c()
typenode1=c()
symbolnode2=c()
typenode2=c()
for (i in 1:length(row.names(allpair))) {
  ID1=as.character(allpair[i,1])
  SYMBOL1=as.character(key[ID1,]$name)
  Type1=as.character(key[ID1,]$type)
  symbolnode1=c(symbolnode1,SYMBOL1)
  typenode1=c(typenode1,Type1)
  
  ID2=as.character(allpair[i,3])
  SYMBOL2=as.character(key[ID2,]$name)
  Type2=as.character(key[ID2,]$type)
  symbolnode2=c(symbolnode2,SYMBOL2)
  typenode2=c(typenode2,Type2)
  }

allpair=data.frame(allpair,symbolnode1,typenode1,symbolnode2,typenode2)
#filter with p<0.05 (or adjusted p<0.05 in the sensitivity analysis)
allpair=allpair[allpair$p<0.05,]
write.table(allpair,file = paste0('~/repo/lncmr_pair_',lyz,'.txt'),sep = '\t',row.names = F)
}




#GO BP
allpair=read.table('~/repo/lncmr_pair_turquoise.txt',sep = '\t',header = T)
#allpair=read.table('~/repo/lncmr_pair_yellow.txt',sep = '\t',header = T)
#allpair=read.table('~/repo/lncmr_pair_magenta.txt',sep = '\t',header = T)
#allpair=read.table('~/repo/lncmr_pair_tan.txt',sep = '\t',header = T)
library(clusterProfiler)
library(org.Hs.eg.db)
ego <- enrichGO(
    gene          = allpair$symbolnode2,
    keyType = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = FALSE)
barplot(ego, showCategory = 20)



