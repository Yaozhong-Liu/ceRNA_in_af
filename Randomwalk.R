#Random Walk with Restart on multiple Networks
library(RandomWalkRestartMH)
library(igraph)

key=read.table("~/repo/key.txt",sep = '\t',header = T)
pair1=read.table('~/repo/lncmr_pair_turquoise.txt',sep = '\t',header = T)
pair2=read.table('~/repo/lncmr_pair_yellow.txt',sep = '\t',header = T)
pair3=read.table('~/repo/lncmr_pair_magenta.txt',sep = '\t',header = T)
pair4=read.table('~/repo/lncmr_pair_tan.txt',sep = '\t',header = T)
allpair=rbind(pair1,pair2,pair3,pair4)
#allpair=allpair[allpair$adjustedp<0.05,] can be used in sensitivity analysis

pro=read.table('G:/生信分析/ceRNA/2020.12/Disease gene database/pro.txt',sep='\t',header = T)
SeedGene=intersect(pro$gene_name,actors[actors$type=='protein_coding',]$name)
#Hypergeometric test
n1=length(actors[actors$type=='protein_coding',]$name)
annotation=read.table("G:/生信分析/ceRNA/Probe_annotation/annotation.txt",sep = '\t',header=T)
n2=length(annotation[annotation$gene_biotype=='protein_coding',]$gene_name)
n3=length(pro$gene_name)
phyper(length(SeedGene)-1,n3,n2-n3,n1,lower.tail = F)


#Random Walk with Restart on two layer Network

#First layer
relations1=allpair[,c('symbolnode1','symbolnode2')]
colnames(relations1)=c('from','to')
actors1=key[key$name%in%c(allpair$symbolnode1,allpair$symbolnode2),]
actors1=actors1[,c('name','type','moduleColors')]
g1 <- graph_from_data_frame(relations1, directed=FALSE, vertices=actors1)

#Second layer
validapair1=data.table::fread('G:/生信分析/ceRNA/2020.12/Pair/cardiac_muscle_top')
summary(validapair1$V3)#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#                           0.1000  0.1001  0.1028  0.1247  0.1177  0.9999  


library(org.Hs.eg.db);
library(clusterProfiler)
dim(actors[actors$type=='lncRNA',])
dim(actors[actors$type=='protein_coding',])

gene.df <- bitr(actors$name, fromType = "SYMBOL", 
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db) 
validepair=validapair1[validapair1$V1%in%gene.df$ENTREZID&validapair1$V2%in%gene.df$ENTREZID,]
validepair$V1=as.character(validepair$V1)
validepair$V2=as.character(validepair$V2)
validepair$ENTREZID=validepair$V1
validepair=merge(validepair,gene.df,by='ENTREZID')
validepair$ENTREZID=validepair$V2
validepair=merge(validepair,gene.df,by='ENTREZID')
validepair=validepair[order(validepair$V3,decreasing = T),]
relations2=validepair[validepair$V3>0.1028,c(5,6)]
colnames(relations2)=c('from','to')
actors2=key[key$name%in%unique(c(relations2$from,relations2$to)),]
actors2=actors2[,c('name','type','moduleColors')]
g2 <- graph_from_data_frame(relations2, directed=FALSE, vertices=actors2)

## We create a 2-layers Multiplex object
Multiplex <- create.multiplex(g1,g2,Layers_Name=c("cerna","valide"))
Multiplex
AdjMatrix_multiplex <- compute.adjacency.matrix(Multiplex)
AdjMatrixNorm_multiplex <- normalize.multiplex.adjacency(AdjMatrix_multiplex)




#set parametyer
ratio=nrow(relations)/nrow(relations2)
c(2/(1+ratio),2*ratio/(1+ratio))
RWR_multiplex_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_multiplex, Multiplex,SeedGene,r=0.7,tau=c(1.78,0.22))#tau=c(1,1)#primary
results=RWR_multiplex_Results$RWRM_Results
colnames(results)=c('name','score')
results=merge(results,actors1,by='name');results=results[order(results$score,decreasing = T),]
write.table(d,file='~/repo/MRWR.txt',sep = '\t')



