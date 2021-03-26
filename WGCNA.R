#read matrix and select the top 5000 genes by variance
rt=read.table('~/repo/rt.txt',sep = '\t',header = T)
rt$var=apply(rt[,1:235],1,var);rt=rt[order(rt$var,decreasing = T),]
datE=rt[1:5000,]
#read phenotype
phenotype=read.table('~/repo/phenotypes_euro.txt',sep = ',',header = T)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF/AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR/SR')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF/SR')

datExpr=datE[,1:235]
library(WGCNA)
datExpr=t(datExpr)
datExpr=data.frame(datExpr)
gsg = goodSamplesGenes(datExpr)
if (!gsg$allOK) {
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

library("doParallel")
## Select the best power parameter
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,blockSize = 20000,networkType="signed",RsquaredCut = 0.9)
sft

#plot scale independence
par(mfrow=c(1,2))
plot(
  sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,
  cex=0.9,
  col="red"
)
abline(h=0.9, col="red")
# plot mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=0.9, col="red")
dev.off()

# run WGCNA
cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power =sft$powerEstimate ,
  maxBlockSize = 20000,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  deepSplit = 2,#2
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs =FALSE,
  verbose = 3)
cor<-stats::cor
table(net$colors)
table(net$unmergedColors)

#visualization
unmergedColors = labels2colors(net$unmergedColors)
mergedColors   = labels2colors(net$colors)
plotDendroAndColors(
  net$dendrograms[[1]],
  cbind(unmergedColors[net$blockGenes[[1]]], mergedColors[net$blockGenes[[1]]]),
  c("Dynamic Tree Cut" , "Merged colors"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

moduleColors = labels2colors(net$colors)
MEs = net$MEs
MEs_col = MEs
library(stringr)
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

#plot eigengene adjacency heatmap
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = F, 
                      xLabelsAngle = 90)


# relationship between clinical information and modules
# binarize it into pairwise indicators
out = binarizeCategoricalVariable(phenotype$AtrialRhythm,
                                  includePairwise = TRUE,
                                  includeLevelVsAll = FALSE);
traitData=data.frame(1-out)
colnames(traitData)=c('AF/AF vs AF/SR persistence','AF/AF vs SR/SR susceptibility and persistence','AF/SR vs SR/SR susceptibility')
traitData=traitData[,c(3,1)]
row.names(traitData)=row.names(datExpr)
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#calculate Pearson's correlation
CORR=data.frame(a=rep(1,length(colnames(MEs))))
PVALUE=data.frame(a=rep(1,length(colnames(MEs))))
for (i in 1:length(colnames(traitData))) {
  R=c()
  P=c()
  for (k in 1:length(colnames(MEs))) {
    r=cor.test(traitData[[i]],MEs[[k]],method = 'pearson')$estimate
    p=cor.test(traitData[[i]],MEs[[k]],method = 'pearson')$p.value
    R=c(R,r)
    P=c(P,p)
  }
  R=as.data.frame(R)
  row.names(R)=colnames(MEs)
  colnames(R)=paste0('GS',colnames(traitData)[[i]])
  CORR=cbind(CORR,R)
  
  P=as.data.frame(P)
  row.names(P)=colnames(MEs)
  colnames(P)=paste0('GS',colnames(traitData)[[i]])
  PVALUE=cbind(PVALUE,P)
}
CORR$a=NULL
PVALUE$a=NULL

#plot
moduleTraitCor=as.matrix(CORR)
moduleTraitPvalue=as.matrix(PVALUE)

textMatrix =  paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1),
  ")",
  sep = ""
)
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs), 
               cex.lab = 0.6, 
               ySymbols = colnames(MEs), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.6, zlim = c(-1,1),
               xLabelsAngle = 0,
               textAdj = c(0.5, 0.5),
               xLabelsAdj = 0.5,
               font.lab.x = 2,
               main = paste("Module-trait relationships"))
dev.off()



# results of WGCNA
#calculate connectivity
ADJ1 = abs(cor(datExpr,use = "p"))^6
Alldegrees = intramodularConnectivity(ADJ1,moduleColors)
#calculate module member-ship
datkME = as.data.frame(cor(datExpr, MEs, use ="p"))
head(datkME)
#write table
connectivity=as.data.frame(Alldegrees);
MM=as.data.frame(datkME);color=as.data.frame(moduleColors)
key=cbind(connectivity,MM,color)
row.names(key)=row.names(datE)
key$type=datE$type
key$name=datE$name
write.table(key,file = paste0('~/repo/key.txt'),sep = '\t')
write.table(MEs,file = paste0('~/repo/MEs.txt'),sep = '\t')


#GO BP 
library(clusterProfiler)
library(org.Hs.eg.db)
name=colnames(MEs)
library(stringr)
name=str_remove(name,'ME')
BP=c()
for (i in name) {
  
  mRNA=key[key$moduleColors==i&key$type=='protein_coding',]
  ego <- enrichGO(
    gene          = mRNA$name,
    keyType = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = FALSE)
  ego=as.data.frame(ego)
  ego$module=i
  ego$number_of_lncRNA=length(key[key$moduleColors==i&key$type=='lncRNA',]$name)
  ego$number_of_proteincoding=length(mRNA$name)
  BP=rbind(BP,ego[1:10,])
}
write.csv(BP,file = '~/repo/GOBP.csv')