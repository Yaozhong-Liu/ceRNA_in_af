rt=read.table('~/repo/rt.txt',sep = '\t',header = T)
lncRNA=read.table('~/repo/MRWR.txt')
lncRNA=lncRNA[lncRNA$type=='lncRNA',]
prepare=rt[rt$name%in%lncRNA$name,]

phenotype=read.table('~/repo/phenotypes_euro.txt',sep = ',',header = T)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF/AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR/SR')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF/SR')

row.names(prepare)=prepare$name
prepare=prepare[,row.names(phenotype)]
lncRNA$AUC_sus='happy'
lncRNA$AUC_per='happy'
lncRNA$AUC_susplusper='happy'
prepare=prepare[lncRNA$name,]


library(pROC)
for (i in 1:dim(prepare)[1]) {
  phenotype1=phenotype[phenotype$AtrialRhythm!='AF/SR',]
  lnc.roc <- roc(response = as.numeric(as.factor(phenotype1$AtrialRhythm)),predictor =as.numeric(prepare[i,row.names(phenotype1)]))
  lncRNA[i,]$AUC_susplusper=lnc.roc$auc
  
  phenotype1=phenotype[phenotype$AtrialRhythm!='SR/SR',]
  lnc.roc <- roc(response = as.numeric(as.factor(phenotype1$AtrialRhythm)),predictor =as.numeric(prepare[i,row.names(phenotype1)]))
  lncRNA[i,]$AUC_per=lnc.roc$auc
  
  phenotype1=phenotype[phenotype$AtrialRhythm!='AF/AF',]
  lnc.roc <- roc(response = as.numeric(as.factor(phenotype1$AtrialRhythm)),predictor =as.numeric(prepare[i,row.names(phenotype1)]))
  lncRNA[i,]$AUC_sus=lnc.roc$auc
  }


write.table(lncRNA,file = '~/repo/lncRNAAUC.txt',sep = '\t')




