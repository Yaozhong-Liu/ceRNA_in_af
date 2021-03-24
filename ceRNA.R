annotation=read.table("G:/生信分析/ceRNA/Probe_annotation/annotation.txt",sep = '\t',header=T)
colnames(annotation)=c('geneId','genesymbol','genebiotype')

#Tarbase8.0
library(dplyr)
Tarbase=read.table('G:/生信分析/lncRNA_miRNA_mRNA database/tarbase/TarBase_v8_download.txt',sep = '\t',header = T)
Tarbase=Tarbase[Tarbase$species=='Homo sapiens',]
Tarbase=Tarbase[,1:3]
Tarbase=distinct(Tarbase)#422893
Tarbase=merge(Tarbase,annotation,by='geneId')#419328
#table(Tarbase$genebiotype)
lnc_Tar=Tarbase[Tarbase$genebiotype=='lncRNA',]#645 pairs
pro_Tar=Tarbase[Tarbase$genebiotype=='protein_coding',]#417860 pairs




#lncbase2.0 experimental validated
lncbase=read.csv('G:/生信分析/lncRNA_miRNA_mRNA database/lncbase/LncBasev2_download.csv',header = F)
library(tidyverse)
lncbase=separate(data = lncbase, col = V1, into = c('geneId','geneName','mirna','species','cell_line','tissue','category','tmethod','positive_negative','direct_indirect','condition'), sep = "\t")
lncbase=lncbase[-1,]
lncbase=lncbase[lncbase$species=='Homo sapiens',]
lncbase=lncbase[,c(1:3)]
lncbase=distinct(lncbase)#62273
lncbase=merge(lncbase,annotation,by='geneId')#47404
#table(lncbase$genebiotype)
lnc_lnc=lncbase[lncbase$genebiotype=='lncRNA',]#25633 pairs
pro_lnc=lncbase[lncbase$genebiotype=='protein_coding',]#906 pairs

Tarbase=rbind(pro_Tar,pro_lnc)
Tarbase=Tarbase[,1:3]
Tarbase=distinct(Tarbase)#418758 pairs

lncbase=rbind(lnc_lnc,lnc_Tar)
lncbase=lncbase[,1:3]
lncbase=distinct(lncbase)#26178 pairs



#co=intersect(unique(lncbase$mirna),unique(Tarbase$mirna))
#lncbase=lncbase[lncbase$mirna%in%co,]
#Tarbase=Tarbase[Tarbase$mirna%in%co,]


#Ready fpr input
key=read.table("G:/生信分析/ceRNA/2020.12/key.txt",sep = '\t',header = T)
MEs=read.table("G:/生信分析/ceRNA/2020.12/MEs.txt",sep = '\t',header = T)
GWAS=read.csv('G:/GWAS/2020.11.20/MY.csv')
AFgene=read.table('G:/生信分析/ceRNA/2020.12/Disease gene database/all.txt',sep='\t',header = T)
AFgene1=names(table(AFgene$x))[table(AFgene$x)>1]
AFgene2=names(table(AFgene$x))[table(AFgene$x)>0]

table(key$moduleColors)
library(stringr)
color=unique(colnames(MEs))
color=str_remove(color,'ME')
color
GO=c()
REA=c()

library('tcltk')
pb <- tkProgressBar("进度","已完成 %", 0, 100)
#u=1:(length(color)-1)
u=c(1,5,16,18)
LYZ=1
for (LYZ in u){
lyz=color[LYZ]
key1=key[key$moduleColors==lyz,]
table(key1$type)
mRNA=key1[key1$type=='protein_coding',]
lncRNA=key1[!key1$type=='protein_coding',]

#lncRNA=lncRNA[lncRNA$kWithin>0.8,]
#mRNA=mRNA[mRNA$kWithin>0.8,]

#Ready for lnc and mrna
library(dplyr)
lnc=lncbase[lncbase$geneId%in%row.names(lncRNA),]
#lnc=distinct(lnc)
mrna=Tarbase[Tarbase$geneId%in%row.names(mRNA),]
#mrna=distinct(mrna)

###################--------------------------------------------------------------noncoding-miRNA-mRNA pairs

cerna=merge(lnc,mrna,by='mirna')
cerna=cerna[,c(2,1,4)]
colnames(cerna)=c('lncID','mirna','mrnaID')
#(unique(cerna$geneName.y))
cerna$mirna=as.character(cerna$mirna)
cerna$lncID=as.character(cerna$lncID)
cerna$mrnaID=as.character(cerna$mrnaID)


if(!dim(cerna)[1]==0){
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

#BH
cerna1=cerna[,c(1,3,4)]
cerna1=distinct(cerna1)
library(qvalue)
cerna1$adjp=p.adjust(cerna1$p, method = 'BH', n = length(cerna1$p))
cerna1$merge=paste(cerna1$lncID,cerna1$mrnaID)
cerna$merge=paste(cerna$lncID,cerna$mrnaID)
final=merge(cerna,cerna1,by='merge')
final=final[,c(2,3,4,8,9)]
colnames(final)=c('lncID','mirna','mrnaID','p','adjusted p')


#####################--------------------------------------------------------------noncoding-noncoding pairs
#library(dplyr)
#noncoding=merge(lnc,lnc,by='mirna')
#noncoding=noncoding[!noncoding$geneId.x==noncoding$geneId.y,]
#noncoding=noncoding[,c(2,1,4)]
#colnames(noncoding)=c('lnc1','mirna','lnc2')


#a=c()
#for (i in 1:length(row.names(noncoding))) {
#  b=as.character(noncoding[i,]$lnc1)
#  c=as.character(noncoding[i,]$lnc2)
#  e=as.character(noncoding[i,]$mirna)
#  d=c(b,c,e)
#  d=sort(d)
#  d=paste(d[1],d[2],d[3])
#  a=c(a,d)
#  }

#noncoding=noncoding[!duplicated(a),]


#hypergeometric test
#noncoding$lnc1=as.character(noncoding$lnc1)
#noncoding$mirna=as.character(noncoding$mirna)
#noncoding$lnc2=as.character(noncoding$lnc2)

#N=length(unique(unique(lncbase$mirna),unique(Tarbase$mirna)))#952
#library(stats)
#pvalue=c()
#for (i in 1:length(row.names(noncoding))) {
#  LNCRNA1=noncoding$lnc1[[i]]
#  LNCRNA2=noncoding$lnc2[[i]]
#  t=noncoding[noncoding$lnc1==LNCRNA1&noncoding$lnc2==LNCRNA2,]
#  I=lncbase[lncbase$geneId==LNCRNA1,];I=unique(as.character(I$mirna))
#  m=lncbase[lncbase$geneId==LNCRNA2,];m=unique(as.character(m$mirna))
#  p=phyper(length(unique(t$mirna))-1,length(m),N-length(m),length(I),lower.tail = F)
#  pvalue=c(pvalue,p)
#}
#noncoding$p=pvalue

#BH
#noncoding1=noncoding[,c(1,3,4)]
#noncoding1=distinct(noncoding1)
#noncoding1$adjp=p.adjust(noncoding1$p, method = 'BH', n = length(noncoding1$p))
#noncoding$merge=paste(noncoding$lnc1,noncoding$lnc2)
#noncoding1$merge=paste(noncoding1$lnc1,noncoding1$lnc2)
#final1=merge(noncoding,noncoding1,by='merge')
#final1=final1[,c(2,3,4,8,9)]
#colnames(final1)=c('lncID1','mirna','lncID2','p','adjusted p')











#####################--------------------------------------------------------------mRNA-mRNA pairs
#mrna=merge(mrna,mrna,by='mirna')
#mrna=mrna[!mrna$geneId.x==mrna$geneId.y,]
#mrna=mrna[,c(2,1,4)]
#colnames(mrna)=c('mrna1','mirna','mrna2')

#a=c()
#for (i in 1:length(row.names(mrna))) {
#  b=as.character(mrna[i,]$mrna1)
#  c=as.character(mrna[i,]$mrna2)
#  e=as.character(mrna[i,]$mirna)
#  d=c(b,c,e)
#  d=sort(d)
#  d=paste(d[1],d[2],d[3])
#  a=c(a,d)
#}

#=mrna[!duplicated(a),]


#hypergeometric test
#mrna$mrna1=as.character(mrna$mrna1)
#mrna$mirna=as.character(mrna$mirna)
#mrna$mrna2=as.character(mrna$mrna2)

#N=length(unique(unique(lncbase$mirna),unique(Tarbase$mirna)))#952
#library(stats)
#pvalue=c()
#for (i in 1:length(row.names(mrna))) {
#  MRNA1=mrna$mrna1[[i]]
#  MRNA2=mrna$mrna2[[i]]
#  t=mrna[mrna$mrna1==MRNA1&mrna$mrna2==MRNA2,]
#  I=Tarbase[Tarbase$geneId==MRNA1,];I=unique(as.character(I$mirna))
#  m=Tarbase[Tarbase$geneId==MRNA2,];m=unique(as.character(m$mirna))
#  p=phyper(length(unique(t$mirna))-1,length(m),N-length(m),length(I),lower.tail = F)
#  pvalue=c(pvalue,p)
#}
#mrna$p=pvalue
#BH
#mrna1=mrna[,c(1,3,4)]
#mrna1=distinct(mrna1)
#mrna1$adjp=p.adjust(mrna1$p, method = 'BH', n = length(mrna1$p))
#mrna$merge=paste(mrna$mrna1,mrna$mrna2)
#mrna1$merge=paste(mrna1$mrna1,mrna1$mrna2)
#final2=merge(mrna,mrna1,by='merge')
#final2=final2[,c(2,3,4,8,9)]
#colnames(final2)=c('mrnaID1','mirna','mrnaID2','p','adjusted p')



#####################--------------------------------------------------------------
final$type=rep('cerna',length(row.names(final)))
write.table(final,file = paste0('G:/生信分析/ceRNA/2020.12/finalpair_',lyz,'.txt'),sep = '\t',row.names = F)

#final=read.table('G:/生信分析/ceRNA/2020.12/finalpair_tan.txt')
#pair=final
#colnames(pair)=c('node1','mirna','node2','p','adjustedp','type')


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
allpair=allpair[allpair$p<0.05,]
#allpair=allpair[,c(7,2,9)]
#colnames(allpair)=c('lncRNA','miRNA','mRNA')
write.csv(allpair,file = 'G:/生信分析/ceRNA/2020.12/triplenetwork_tan.csv',row.names = F)

write.table(allpair,file = paste0('G:/生信分析/ceRNA/2020.12/allpair_',lyz,'.txt'),sep = '\t',row.names = F)
#unique(allpair$symbolnode1);length(unique(allpair$symbolnode1))
#unique(allpair$symbolnode2);length(unique(allpair$symbolnode2))

allpair=read.table('G:/生信分析/ceRNA/2020.12/allpair_turquoise.txt',sep = '\t',header = T)
allpair=read.table('G:/生信分析/ceRNA/2020.12/allpair_yellow.txt',sep = '\t',header = T)
allpair=read.table('G:/生信分析/ceRNA/2020.12/allpair_magenta.txt',sep = '\t',header = T)
allpair=read.table('G:/生信分析/ceRNA/2020.12/allpair_tan.txt',sep = '\t',header = T)
write1=0
#write1=dim(allpair)[1]
if (!write1==0)
{
write2=length(unique(allpair$symbolnode1))
write3=length(unique(allpair$symbolnode2))
a=sort(table(allpair$symbolnode1),decreasing = T)
a=list(a)
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
ego=as.data.frame(ego)
ego$module=lyz
ego$numberofcerna=paste0(write1,' (',write2,' lncRNAs and ',write3,' mRNAs)')
ego$lncrnasortbydegree=paste(a)
ego$originalnumberofrna=paste0(dim(lncRNA)[1],' ',dim(mRNA)[1])
ego$gwas=paste(GWAS[GWAS$gene_id%in%unique(c(as.character(allpair$node1),as.character(allpair$node2))),2],collapse=',')
ego$afgene1=paste(intersect(AFgene1,unique(c(as.character(allpair$symbolnode1),as.character(allpair$symbolnode2)))),collapse = ',')
ego$afgene2=paste(intersect(AFgene2,unique(c(as.character(allpair$symbolnode1),as.character(allpair$symbolnode2)))),collapse = ',')
GO=rbind(GO,ego)




gene.df <- bitr(allpair$symbolnode2, fromType = "SYMBOL", 
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db) 
library(ReactomePA);library(reactome.db)
rea=enrichPathway(
  gene=gene.df$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)
barplot(rea, showCategory = 20)
rea=as.data.frame(rea)
rea$module=lyz
rea$numberofcerna=paste0(write1,' (',write2,' lncRNAs and ',write3,' mRNAs)')
rea$lncrnasortbydegree=paste(a)
rea$originalnumberofrna=paste0(dim(lncRNA)[1],' ',dim(mRNA)[1])
rea$gwas=paste(GWAS[GWAS$gene_id%in%unique(c(as.character(allpair$node1),as.character(allpair$node2))),2],collapse=',')
rea$afgene1=paste(intersect(AFgene1,unique(c(as.character(allpair$symbolnode1),as.character(allpair$symbolnode2)))),collapse = ',')
rea$afgene2=paste(intersect(AFgene2,unique(c(as.character(allpair$symbolnode1),as.character(allpair$symbolnode2)))),collapse = ',')
REA=rbind(REA,rea)
}
}
info=sprintf("已完成 %d%%", round(LYZ*100/length(u)))  
setTkProgressBar(pb, LYZ*100/length(u), sprintf("进度 (%s)", info),info)
}

write.csv(GO,file = 'G:/生信分析/ceRNA/2020.12/WGCNA_GO.csv') 
write.csv(REA,file = 'G:/生信分析/ceRNA/2020.12/WGCNA_REA.csv') 

library(stringr)
color=unique(colnames(MEs))
color=str_remove(color,'ME')

GO=read.csv('G:/生信分析/ceRNA/2020.12/WGCNA_GO.csv',header = T,)
GO=GO[GO$p.adjust<0.05,]
head(GO[GO$module==color[1],])

REA=read.csv('G:/生信分析/ceRNA/2020.12/WGCNA_REA.csv',header = T,)
REA=REA[REA$p.adjust<0.05,]
head(REA[REA$module==color[4],])

lyz=REA[REA$module==color[17],]
#####################--------------------------------------------------------------
#key[Reduce(intersect,list(row.names(key),my$gene_id,lncbase$geneId)),]
Trans=read.csv("G:/生信分析/Integrative omics data/BMC 10.13/REM.csv")
key[key$moduleColors==color[2]&key$type=='lncRNA',]
Trans[Trans$FDR<0.05&Trans$X%in%key[key$moduleColors==color[7]&key$type=='lncRNA',]$name,]




