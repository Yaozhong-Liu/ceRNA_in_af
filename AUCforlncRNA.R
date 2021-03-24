rt=read.table('G:/生信分析/ceRNA/2020.12/rt.txt',sep = '\t',header = T,row.names = 1)
lncRNA=read.table('G:/生信分析/ceRNA/2020.12/MRWR_sensitivity.txt')
lncRNA=lncRNA[lncRNA$type=='lncRNA',]
prepare=rt[rt$name%in%lncRNA$name,]
phenotype=read.table('G:/生信分析/ceRNA/2020.12/GSE69890_RAW/phenotypes_euro.txt',sep = ',',header = T,row.names = 1)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
#phenotype=phenotype[phenotype$AtrialRhythm!='Yes AF/Sinus Rhythm',]
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

#library(pROC)
#for (i in 1:dim(prepare)[1]) {
#  phenotype1=phenotype[phenotype$AtrialRhythm!='AF/SR',]
#  lnc.roc <- roc(response = as.numeric(as.factor(phenotype1$AtrialRhythm)),predictor =as.numeric(prepare[i,row.names(phenotype1)]))
#  lncRNA[i,]$AUC=lnc.roc$auc
#  if(lncRNA[i,]$moduleColors%in%c('turquoise','yellow')){
#    phenotype1=phenotype[phenotype$AtrialRhythm!='SR/SR',]
#    lnc.roc <- roc(response = as.numeric(as.factor(phenotype1$AtrialRhythm)),predictor =as.numeric(prepare[i,row.names(phenotype1)]))
#    lncRNA[i,]$AUCsp=lnc.roc$auc
#  }else{
#    phenotype1=phenotype[phenotype$AtrialRhythm!='AF/AF',]
#    lnc.roc <- roc(response = as.numeric(as.factor(phenotype1$AtrialRhythm)),predictor =as.numeric(prepare[i,row.names(phenotype1)]))
#    lncRNA[i,]$AUCsp=lnc.roc$auc
#  }
#}


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



write.table(lncRNA,file = 'G:/生信分析/ceRNA/2020.12/lncRNAAUC_sensitivity.txt',sep = '\t')




