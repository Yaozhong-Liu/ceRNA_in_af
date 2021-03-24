#BiocManager::install('RandomWalkRestartMH')
##Random Walk with Restart on a Monoplex Network
library(RandomWalkRestartMH)
library(igraph)

#Trans=read.csv("G:/生信分析/Integrative omics data/BMC 10.13/REM.csv")
Trans=read.csv("G:/生信分析/Integrative omics data/GeneMeta/REM.csv")
key=read.table("G:/生信分析/ceRNA/2020.12/key.txt",sep = '\t',header = T)


###cerna
#allpair1=read.table('G:/生信分析/ceRNA/2020.12/predict/allpair_turquoise.txt',sep = '\t',header = T)
allpair1=read.table('G:/生信分析/ceRNA/2020.12/allpair_turquoise.txt',sep = '\t',header = T)
#allpair1=allpair1[allpair1$symbolnode1%in%key[key$MEblue>0.5,]$name&allpair1$symbolnode2%in%key[key$MEblue>0.5,]$name,]
#allpair1=allpair1[allpair1$symbolnode1%in%key[key$adj.P.Val.x<0.05,]$name&allpair1$symbolnode2%in%key[key$adj.P.Val.x<0.05,]$name,]

#allpair2=read.table('G:/生信分析/ceRNA/2020.12/predict/allpair_yellow.txt',sep = '\t',header = T)
allpair2=read.table('G:/生信分析/ceRNA/2020.12/allpair_yellow.txt',sep = '\t',header = T)
#allpair2=allpair2[allpair2$symbolnode1%in%key[key$MEgreen>0.5,]$name&allpair2$symbolnode2%in%key[key$MEgreen>0.5,]$name,]
#allpair2=allpair2[allpair2$symbolnode1%in%key[key$adj.P.Val.x<0.05,]$name&allpair2$symbolnode2%in%key[key$adj.P.Val.x<0.05,]$name,]

#allpair3=read.table('G:/生信分析/ceRNA/2020.12/predict/allpair_magenta.txt',sep = '\t',header = T)
allpair3=read.table('G:/生信分析/ceRNA/2020.12/allpair_magenta.txt',sep = '\t',header = T)
#allpair3=allpair3[allpair3$symbolnode1%in%key[key$MEmagenta>0.5,]$name&allpair3$symbolnode2%in%key[key$MEmagenta>0.5,]$name,]
#allpair3=allpair3[allpair3$symbolnode1%in%key[key$adj.P.Val.y<0.05,]$name&allpair3$symbolnode2%in%key[key$adj.P.Val.y<0.05,]$name,]

#allpair4=read.table('G:/生信分析/ceRNA/2020.12/predict/allpair_tan.txt',sep = '\t',header = T)
allpair4=read.table('G:/生信分析/ceRNA/2020.12/allpair_tan.txt',sep = '\t',header = T)
#allpair4=allpair4[allpair4$symbolnode1%in%key[key$MEtan>0.5,]$name&allpair4$symbolnode2%in%key[key$MEtan>0.5,]$name,]
#allpair4=allpair4[allpair4$symbolnode1%in%key[key$adj.P.Val.y<0.05,]$name&allpair4$symbolnode2%in%key[key$adj.P.Val.y<0.05,]$name,]

#allpair5=read.table('G:/生信分析/ceRNA/2020.12/allpair_green.txt',sep = '\t',header = T)
#allpair5=allpair5[allpair5$symbolnode1%in%key[key$MEgreen>0.5,]$name&allpair5$symbolnode2%in%key[key$MEgreen>0.5,]$name,]
#allpair5=allpair5[allpair5$symbolnode1%in%key[key$adj.P.Val.y<0.05,]$name&allpair5$symbolnode2%in%key[key$adj.P.Val.y<0.05,]$name,]


#allpair6=read.table('G:/生信分析/ceRNA/2020.12/allpair_red.txt',sep = '\t',header = T)
#allpair6=allpair6[allpair6$symbolnode1%in%key[key$MEred>0.5,]$name&allpair6$symbolnode2%in%key[key$MEred>0.5,]$name,]
#allpair6=allpair6[allpair6$symbolnode1%in%key[key$adj.P.Val.x<0.05,]$name&allpair6$symbolnode2%in%key[key$adj.P.Val.x<0.05,]$name,]


#allpair7=read.table('G:/生信分析/ceRNA/2020.12/allpair_greenyellow.txt',sep = '\t',header = T)
#allpair7=allpair7[allpair7$symbolnode1%in%key[key$MEgreenyellow>0.5,]$name&allpair7$symbolnode2%in%key[key$MEgreenyellow>0.5,]$name,]
#allpair7=allpair7[allpair7$symbolnode1%in%key[key$adj.P.Val.y<0.05,]$name&allpair7$symbolnode2%in%key[key$adj.P.Val.y<0.05,]$name,]


#allpair8=read.table('G:/生信分析/ceRNA/2020.12/allpair_tan.txt',sep = '\t',header = T)
#allpair8=allpair8[allpair8$symbolnode1%in%key[key$MEtan>0.5,]$name&allpair8$symbolnode2%in%key[key$MEtan>0.5,]$name,]
#allpair8=allpair8[allpair8$symbolnode1%in%key[key$adj.P.Val.x<0.05,]$name&allpair8$symbolnode2%in%key[key$adj.P.Val.x<0.05,]$name,]

#allpair=rbind(allpair1,allpair2,allpair3,allpair4,allpair5,allpair6,allpair7,allpair8)
allpair=rbind(allpair1,allpair2,allpair3,allpair4)
#allpair=rbind(allpair1)



#allpair=allpair[allpair$adjustedp<0.05,]
relations=allpair[,c(6,8)]
colnames(relations)=c('from','to')
actors=key[key$name%in%c(allpair$symbolnode1,allpair$symbolnode2),]
actors=actors[,c('name','type','moduleColors')]
g <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
cerna=create.multiplex(g,Layers_Name = 'cerna')
cerna
AdjMatrix_cerna=compute.adjacency.matrix(cerna)
AdjMatrixNorm_cerna=normalize.multiplex.adjacency(AdjMatrix_cerna)
#Seed gene
#my=read.csv('G:/GWAS/2020.11.20/MY.csv')
#intersect(my$gene_name,actors$name)
#lyz=intersect(my$gene_name,SeedGene)
#lyz

#my=read.table('G:/生信分析/ceRNA/2020.12/Disease gene database/all.txt',sep='\t',header = T)
pro=read.table('G:/生信分析/ceRNA/2020.12/Disease gene database/pro.txt',sep='\t',header = T)


#my=names(table(my$x))[table(my$x)>0]
#actors[actors$name%in%my&actors$type=='lncRNA',]
#SeedGene=intersect(my,actors[actors$type=='protein_coding',]$name)
SeedGene=intersect(pro$gene_name,actors[actors$type=='protein_coding',]$name)

###超几何分布
#n1=length(actors[actors$type=='protein_coding',]$name)
#annotation=read.table("G:/生信分析/ceRNA/Probe_annotation/annotation.txt",sep = '\t',header=T)
#n2=length(annotation[annotation$gene_biotype=='protein_coding',]$gene_name)
#n3=length(pro$gene_name)
#phyper(length(SeedGene)-1,n3,n2-n3,n1,lower.tail = F)
#actors[actors$name%in%SeedGene,]
#relations[relations$to%in%SeedGene,]


#RWR_cerna_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_cerna,cerna,SeedGene,r=0.7)
#RWR_cerna_Results
#a=RWR_cerna_Results$RWRM_Results
#colnames(a)=c('name','score')
#b=merge(a,actors,by='name');b=b[order(b$score,decreasing = T),]
#b[1:10,]
#Trans[Trans$X%in%b[1:5,]$name,]


###Plot
#TopResults <-create.multiplexNetwork.topResults(RWR_cerna_Results,cerna,k=3)
#par(mar=c(0.1,0.1,0.1,0.1))
#plot(TopResults, vertex.label.color="black",vertex.frame.color="#ffffff",vertex.size= 5, 
#     edge.curved=.2,vertex.color = ifelse(igraph::V(TopResults)$name%in%SeedGene,"yellow","#00CCFF"), edge.color="blue",edge.width=0.8)




###########
validapair1=data.table::fread('G:/生信分析/ceRNA/2020.12/Pair/cardiac_muscle_top')
#validapair1=read.table('G:/生信分析/ceRNA/2020.12/Pair/heart_top')
summary(validapair1$V3)#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#                           0.1000  0.1001  0.1028  0.1247  0.1177  0.9999  
#validapair1=validapair1[validapair1$V3>0.1065,]


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


#validepair[validepair$SYMBOL.x=='C15orf56'|validepair$SYMBOL.y=='C15orf56',]
#validepair[validepair$SYMBOL.x=='LINC00964'|validepair$SYMBOL.y=='LINC00964',]
#validepair[validepair$SYMBOL.x=='MIAT'|validepair$SYMBOL.y=='MIAT',]


relations2=validepair[validepair$V3>0.1028,c(5,6)]
#relations2=validepair[1:round(nrow(validepair)/2),c(5,6)]
colnames(relations2)=c('from','to')
actors2=key[key$name%in%unique(c(relations2$from,relations2$to)),]
actors2=actors2[,c('name','type','moduleColors')]


g1 <- graph_from_data_frame(relations2, directed=FALSE, vertices=actors2)
## We create a 2-layers Multiplex object
Multiplex <- create.multiplex(g,g1,Layers_Name=c("cerna","valide"))
Multiplex
AdjMatrix_multiplex <- compute.adjacency.matrix(Multiplex)
AdjMatrixNorm_multiplex <- normalize.multiplex.adjacency(AdjMatrix_multiplex)


ratio=nrow(relations)/nrow(relations2)
c(2/(1+ratio),2*ratio/(1+ratio))
#?Random.Walk.Restart.Multiplex
#RWR_multiplex_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_multiplex, Multiplex,SeedGene,r=0.7,tau=c(1.8,0.2))#tau=c(1,1)#sensitivity
RWR_multiplex_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_multiplex, Multiplex,SeedGene,r=0.7,tau=c(1.78,0.22))#tau=c(1,1)#primary
c=RWR_multiplex_Results$RWRM_Results
colnames(c)=c('name','score')
d=merge(c,actors,by='name');d=d[order(d$score,decreasing = T),]
key[key$name%in%d[1:3,]$name,]
d[1:10,]
dim(d[d$type=='lncRNA',])
write.table(d,file='G:/生信分析/ceRNA/2020.12/MRWR_sensitivity.txt',sep = '\t')





#relations[relations$from=='LINC00964',]
Trans[Trans$X%in%d[1:3,]$name,]


###Plot
TopResults <-create.multiplexNetwork.topResults(RWR_multiplex_Results,Multiplex,k=2)
par(mar=c(0.1,0.1,0.1,0.1))

plot(TopResults, vertex.label.color="black",vertex.frame.color="#ffffff", vertex.size= 5,
     edge.curved= ifelse(E(TopResults)$type == "valide", 0.4,0),
     vertex.color = ifelse(igraph::V(TopResults)$name %in%SeedGene,"yellow","#00CCFF"),edge.width=0.8,
     edge.color=ifelse(E(TopResults)$type == "valide","blue","red")
     )



###Top 50%
candidate=actors[actors$type=='lncRNA',]$name
n=length(candidate)


#MM
MM=key[key$name%in%candidate,]
MM=MM[order(MM$kWithin,decreasing = T),c(2,29)]

#RWRMscore
RWRMscore=d[d$name%in%candidate,]

#
degree=names(table(allpair$symbolnode1))[order(table(allpair$symbolnode1),decreasing = T)]


RWRscore=b[b$name%in%candidate,]





A=intersect(MM[1:10,]$name,RWRMscore[1:10,]$name)
B=intersect(A,RWRscore[1:10,]$name)
#B=intersect(B,degree[1:10])
B



Trans[Trans$X=='LINC00964',]




mRNA=allpair[allpair$symbolnode1=='LINC00964',]
library(clusterProfiler)
#BiocManager::install('fgsea')
library(org.Hs.eg.db)

ego <- enrichGO(
                gene = mRNA$symbolnode2,
                keyType = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)
barplot(ego, showCategory = 10)
head(as.data.frame(ego))




?enrichGO



