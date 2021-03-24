setwd('G:/生信分析/ceRNA/2020.12')
datExpr=datE[,1:235]
#boxplot(datExpr,outline = F)
library(WGCNA)
datExpr=t(datExpr)
datExpr=data.frame(datExpr)
gsg = goodSamplesGenes(datExpr)
if (!gsg$allOK) {
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#D=dist(datExpr,diag = FALSE,upper = FALSE)
#hc=hclust(D)
#plot(hc,phenotype$AtrialRhythm)
#rect.hclust(hc,2)
###


library("doParallel")
## Select the best power parameter 0.85 or 0.9
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
sft$powerEstimate#12
dev.off()

## RUN WGCNA(Signed or unsigned) and save Toms or not
cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power =sft$powerEstimate ,#sft$powerEstimate
  maxBlockSize = 20000,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,#30
  reassignThreshold = 0,
  mergeCutHeight = 0.25,#0.25
  deepSplit = 2,#2
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs =FALSE,
  #saveTOMFileBase = "AF-blockwise",
  verbose = 3)
cor<-stats::cor
table(net$colors)
table(net$unmergedColors)






###select parameter
  validapair1=read.table('G:/生信分析/ceRNA/2020.12/Pair/heart_top')
  validapair1$pair=paste(validapair1$V1,validapair1$V2)
  #validapair2=read.table('G:/生信分析/ceRNA/2020.12/Pair/cardiac_muscle_top')
  #validapair2$pair=paste(validapair2$V1,validapair2$V2)
  Parameter=read.table('G:/生信分析/ceRNA/2020.12/Parameter set/parameter.txt',sep = '\t',header = T)
  nomodule=c()
  noingrey=c()
  hsmsp1=c()
  hsmnsp1=c()
  hnsmsp1=c()
  hnsmnsp1=c()
  hscorewithgrey=c()
  hsmsp2=c()
  hsmnsp2=c()
  hnsmsp2=c()
  hnsmnsp2=c()
  hscorewithoutgrey=c()
  #csmsp1=c()
  #csmnsp1=c()
  #cnsmsp1=c()
  #cnsmnsp1=c()
  #cscorewithgrey=c()
  #csmsp2=c()
  #csmnsp2=c()
  #cnsmsp2=c()
  #cnsmnsp2=c()
  #cscorewithoutgrey=c()
  library('tcltk')
  pb <- tkProgressBar("进度","已完成 %", 0, 100)
  u=1:dim(Parameter)[1]
  for (LYZ in u) {
    pow=Parameter[LYZ,1]
    size=Parameter[LYZ,2]
    split=Parameter[LYZ,3]
    cuttree=Parameter[LYZ,4]
  
  
    ## RUN WGCNA(Signed or unsigned) and save Toms or not
    cor <- WGCNA::cor
    net = blockwiseModules(
      datExpr,
      power =pow ,#sft$powerEstimate
      maxBlockSize = 20000,
      TOMType = "signed",
      networkType = "signed",
      minModuleSize = size,#30
      reassignThreshold = 0,
      mergeCutHeight = cuttree,#0.25
      deepSplit = split,#2
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      saveTOMs =FALSE,
      #saveTOMFileBase = "AF-blockwise",
      verbose = 3)
    cor<-stats::cor
    table(net$colors)
    table(net$unmergedColors)
    
    mergedColors= labels2colors(net$colors)
    
    #模块数，灰色模块数
    nomodule=c(nomodule,length (unique (mergedColors)))#include grey
    noingrey=c(noingrey,length (mergedColors[mergedColors=='grey']))
    #模块数，灰色模块数
    
   
    #BiocManager::install("org.Hs.eg.db")
    library(org.Hs.eg.db);
    library(clusterProfiler)
    
    gene.df <- bitr(row.names(datE), fromType = "ENSEMBL", 
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db) 
    
  
    ######################################################################################validapair1
    copair=datE
    #copair$modulecolor=mergedColors
    copair$modulecolor=copair$moduleColors
    copair$ENSEMBL=rownames(copair)
    copair=merge(copair,gene.df,by='ENSEMBL')
    
    copairall=copair[copair$ENTREZID%in%unique(c(validapair1$V1,validapair1$V2)),]
    copairall=copairall[!duplicated(copairall$ENSEMBL),]
    copairall=copairall[order(copairall$ENTREZID),]
    
    
    copair=copair[!copair$modulecolor=='grey',]
    copair=copair[copair$ENTREZID%in%unique(c(validapair1$V1,validapair1$V2)),]
    copair=copair[!duplicated(copair$ENSEMBL),]
    copair=copair[order(copair$ENTREZID),]
    
    Modules <- unique (copair$modulecolor)
    module_pairs <- c ()
    for (module in 1 : length (Modules)) { 
      subset <- copair[copair$modulecolor==Modules[module], ]
      Pairs <- as.data.frame (t (combn (subset$ENTREZID, 2)))
      Pairs$module <- Modules [ module ] 
      module_pairs [[ module ]] <- Pairs
    }
    
    #module pair not grey
    module_pairs <- as.data.frame (do.call (rbind, module_pairs))
    head (module_pairs)
    dim (module_pairs)
    
    #all pair (including grey genes)
    All_pairs=as.data.frame (t (combn (copairall$ENTREZID, 2)))
    dim(All_pairs)
    
    #all pair (without grey genes)
    All_pairs_withoutgrey=as.data.frame (t (combn (copair$ENTREZID, 2)))
    dim(All_pairs_withoutgrey)
    
    #set pair
    module_pairs$pair=paste(module_pairs$V1,module_pairs$V2)
    All_pairs$pair=paste(All_pairs$V1,All_pairs$V2)
    All_pairs_withoutgrey$pair=paste(All_pairs_withoutgrey$V1,All_pairs_withoutgrey$V2)
    
    #calculate and write
    SM_SP1=dim(module_pairs[module_pairs$pair%in%validapair1$pair,])[1]
    SM_nSP1=dim(module_pairs)[1]-SM_SP1
    nSM_SP1=dim(All_pairs[All_pairs$pair%in%validapair1$pair,])[1]-SM_SP1
    nSM_nSP1=dim(All_pairs)[1]-nSM_SP1-SM_SP1-SM_nSP1
    hsmsp1=c(hsmsp1,SM_SP1)
    hsmnsp1=c(hsmnsp1,SM_nSP1)
    hnsmsp1=c(hnsmsp1,nSM_SP1)
    hnsmnsp1=c(hnsmnsp1,nSM_nSP1)
    hscorewithgrey=c(hscorewithgrey,(SM_SP1/nSM_SP1)*(nSM_nSP1/SM_nSP1))
    
    
    
    SM_SP2=SM_SP1
    SM_nSP2=SM_nSP1
    nSM_SP2=dim(All_pairs_withoutgrey[All_pairs_withoutgrey$pair%in%validapair1$pair,])[1]-SM_SP2
    nSM_nSP2=dim(All_pairs_withoutgrey)[1]-nSM_SP2-SM_SP2-SM_nSP2
    hsmsp2=c(hsmsp2,SM_SP2)
    hsmnsp2=c(hsmnsp2,SM_nSP2)
    hnsmsp2=c(hnsmsp2,nSM_SP2)
    hnsmnsp2=c(hnsmnsp2,nSM_nSP2)
    hscorewithoutgrey=c(hscorewithoutgrey,(SM_SP2/nSM_SP2)*(nSM_nSP2/SM_nSP2))
    ######################################################################################validapair1
    
    
    info=sprintf("已完成 %d%%", round(i*100/length(u)))  
    setTkProgressBar(pb, i*100/length(u), sprintf("进度 (%s)", info),info)
    }
  output=data.frame(Parameter,nomodule,noingrey,hsmsp1,hsmnsp1,hnsmsp1,hnsmnsp1,hscorewithgrey,hsmsp2,hsmnsp2,hnsmsp2,hnsmnsp2,hscorewithoutgrey)
  write.csv(output,file = 'output_500_4500.csv')
  output=read.csv('output.csv')
  
  
  
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
  
  
  ## Relationship between clinical information and modules
  
  #library(WGCNA)
  # Binarize it into pairwise indicators
  out = binarizeCategoricalVariable(phenotype$AtrialRhythm,
                                    includePairwise = TRUE,
                                    includeLevelVsAll = FALSE);
  # Print the variable and the indicators
  data.frame(phenotype$AtrialRhythm,out);
  
  
  #####
  #AF.SR.AF.AF=phenotype[!phenotype$AtrialRhythm=='SR/SR',]
  #model1=model.matrix(~sexFemale+AtrialRhythm,data = AF.SR.AF.AF)
  #MEs1=MEs[row.names(AF.SR.AF.AF),]
  #  R=c()
  #  P=c()
    
  #  for (i in 1:length(colnames(MEs1))) {
  #    fit=lm(MEs1[[i]]~model1)
  #    r=fit$coefficients[4]
  #    p=coef(summary(fit))[3,4] 
  #   R=c(R,r)
  #    P=c(P,p)
  #       }
  #  R=as.data.frame(R)
  #  row.names(R)=colnames(MEs)
  #  P=as.data.frame(P)
  #  row.names(P)=colnames(MEs)
  ####
  
  
  
  
  #traitData=data.frame(4-as.numeric(as.factor(phenotype$AtrialRhythm)),1-out)
  traitData=data.frame(1-out)
  colnames(traitData)=c('AF/AF vs AF/SR persistence','AF/AF vs SR/SR susceptibility and persistence','AF/SR vs SR/SR susceptibility')
  traitData=traitData[,c(3,1)]
  traitData
  
  #colnames(traitData)=c('Progression','AF/AF vs AF/SR persistence','AF/AF vs SR/SR susceptibility and persistence','AF/SR vs SR/SR susceptibility')
  row.names(traitData)=row.names(datExpr)
  corType = "pearson"
  corFnc = ifelse(corType=="pearson", cor, bicor)
  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  robustY = ifelse(corType=="pearson",T,F)
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  
  
  ################################################################################
  traitData=data.frame(traitData)
  MEs=MEs
  
  #gene signigicancer
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
  
  ########################################################################################################
  
  #moduleTraitCor = cor(
  #  MEs,
  #  traitData,
  #  use = "p"
  #)
  
  
  
  ##moduleTraitPvalue = corPvalueStudent(
  #  moduleTraitCor,
  #  nSamples
  #)
  
  
  
  ########################################################################################################
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
  #plot
  
  #pdf(file=paste('wgcna',LYZ,'.pdf'),width=6, height=10)
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
  
  
  datE$moduleColors=mergedColors
  key=datE[,236:239]
  table(key$moduleColors)
  key[key$name=='MIAT',]
  lyz=key[key$moduleColors=='tan',]
  ## Results of WGCNA
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
  
  
  allDiff1=read.table('G:/生信分析/ceRNA/2020.12/alldiff_persistence.txt',sep = '\t',header = T)
  allDiff2=read.table('G:/生信分析/ceRNA/2020.12/alldiff_susceptible.txt',sep = '\t',header = T)
  
  key$ID=rownames(key)
  allDiff1$ID=row.names(allDiff1)
  key=merge(key,allDiff1,by='ID')
  allDiff2$ID=row.names(allDiff2)
  key=merge(key,allDiff2,by='ID')
  ##xpersistence ysusceptibility
  row.names(key)=key$ID
  key$ID=NULL
  
  lyz=key[key$moduleColors=='yellow',]
  key[key$name=='MIAT',]
  #write.table(key,file = paste0('G:/生信分析/ceRNA/2020.12/key',LYZ,'.txt'),sep = '\t')
  #write.table(MEs,file = paste0('G:/生信分析/ceRNA/2020.12/MEs',LYZ,'.txt'),sep = '\t')
  
  write.table(key,file = paste0('G:/生信分析/ceRNA/2020.12/key.txt'),sep = '\t')
  write.table(MEs,file = paste0('G:/生信分析/ceRNA/2020.12/MEs.txt'),sep = '\t')
  
  
  #info=sprintf("已完成 %d%%", round(LYZ*100/length(u)))  
  #setTkProgressBar(pb, LYZ*100/length(u), sprintf("进度 (%s)", info),info)
  #}


table(key$moduleColors)
my=read.table('G:/生信分析/ceRNA/2020.12/Disease gene database/all.txt',sep='\t',header = T)
key$afgene=as.character(key$name%in%my$x)

########################################################################################################
key=read.table("G:/生信分析/ceRNA/2020.12/key.txt",sep = '\t',header = T)

index=unique(key$moduleColors)
for (i in 1:length(index)) {
  print(index[i])
  print(table(key[key$moduleColors==index[i],paste0('ME',index[i])]>=0))
  
}



#my=read.csv('G:/GWAS/2020.11.20/MY.csv')
#b=key[row.names(key)%in%my$gene_id,]
#b=key[key$name%in%my1,]


Trans=read.csv("G:/生信分析/Integrative omics data/BMC 10.13/REM.csv")
Trans=Trans[Trans$FDR<0.05,]
#table(a$moduleColors)
Trans[Trans$X%in%key[key$moduleColors=='blue',]$name,c(1,18,19)]


#a=read.table('G:/生信分析/ceRNA/2020.12/proteomic.txt',sep = '\t',header = T)
#a=a[a$P.value<0.05,]
#a=a[a$Gene.name%in%key[key$moduleColors=='turquoise',]$name,]
#c=key[key$name%in%a$Gene.name,]

key=read.table("G:/生信分析/ceRNA/2020.12/key.txt",sep = '\t',header = T)
MEs=read.table("G:/生信分析/ceRNA/2020.12/MEs.txt",sep = '\t',header = T)
#-------------------------------------------------GO富集分析
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
write.csv(BP,file = 'G:/生信分析/ceRNA/2020.12/BP.csv')



######################################################################################validapair2
copair=datE
copair$modulecolor=mergedColors

#模块数，灰色模块数
nomodule=c(nomodule,length (unique (copair$modulecolor)))#include grey
noingrey=c(noingrey,dim(copair[copair$modulecolor=='grey',])[1])
#模块数，灰色模块数

copair$ENSEMBL=rownames(copair)
copair=merge(copair,gene.df,by='ENSEMBL')

copairall=copair[copair$ENTREZID%in%unique(c(validapair2$V1,validapair2$V2)),]
copairall=copairall[!duplicated(copairall$ENSEMBL),]
copairall=copairall[order(copairall$ENTREZID),]


copair=copair[!copair$modulecolor=='grey',]
copair=copair[copair$ENTREZID%in%unique(c(validapair2$V1,validapair2$V2)),]
copair=copair[!duplicated(copair$ENSEMBL),]
copair=copair[order(copair$ENTREZID),]

Modules <- unique (copair$modulecolor)
module_pairs <- c ()
for (module in 1 : length (Modules)) { 
  subset <- copair[copair$modulecolor==Modules[module], ]
  Pairs <- as.data.frame (t (combn (subset$ENTREZID, 2)))
  Pairs$module <- Modules [ module ] 
  module_pairs [[ module ]] <- Pairs
}

#module pair not grey
module_pairs <- as.data.frame (do.call (rbind, module_pairs))
head (module_pairs)
dim (module_pairs)

#all pair (including grey genes)
All_pairs=as.data.frame (t (combn (copairall$ENTREZID, 2)))
dim(All_pairs)

#all pair (without grey genes)
All_pairs_withoutgrey=as.data.frame (t (combn (copair$ENTREZID, 2)))
dim(All_pairs_withoutgrey)

#set pair
module_pairs$pair=paste(module_pairs$V1,module_pairs$V2)
All_pairs$pair=paste(All_pairs$V1,All_pairs$V2)
All_pairs_withoutgrey$pair=paste(All_pairs_withoutgrey$V1,All_pairs_withoutgrey$V2)

#calculate and write
SM_SP1=dim(module_pairs[module_pairs$pair%in%validapair2$pair,])[1]
SM_nSP1=dim(module_pairs)[1]-SM_SP1
nSM_SP1=dim(All_pairs[All_pairs$pair%in%validapair2$pair,])[1]-SM_SP1
nSM_nSP1=dim(All_pairs)[1]-nSM_SP1-SM_SP1-SM_nSP1
csmsp1=c(csmsp1,SM_SP1)
csmnsp1=c(csmnsp1,SM_nSP1)
cnsmsp1=c(cnsmsp1,nSM_SP1)
cnsmnsp1=c(cnsmnsp1,nSM_nSP1)
cscorewithgrey=c(cscorewithgrey,(SM_SP1/nSM_SP1)*(nSM_nSP1/SM_nSP1))



SM_SP2=SM_SP1
SM_nSP2=SM_nSP1
nSM_SP2=dim(All_pairs_withoutgrey[All_pairs_withoutgrey$pair%in%validapair1$pair,])[1]-SM_SP2
nSM_nSP2=dim(All_pairs_withoutgrey)[1]-nSM_SP2-SM_SP2-SM_nSP2
csmsp2=c(csmsp2,SM_SP2)
csmnsp2=c(csmnsp2,SM_nSP2)
cnsmsp2=c(cnsmsp2,nSM_SP2)
cnsmnsp2=c(cnsmnsp2,nSM_nSP2)
cscorewithoutgrey=c(cscorewithoutgrey,(SM_SP2/nSM_SP2)*(nSM_nSP2/SM_nSP2))
######################################################################################validapair2

