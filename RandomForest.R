allpair1=read.table('~/repo/lncmr_turquoise.txt',sep = '\t',header = T)
allpair1=allpair1[allpair1$symbolnode1%in%c('LINC00964'),]

#allpair2=read.table('G:/ÉúÐÅ·ÖÎö/ceRNA/2020.12/allpair_tan.txt',sep = '\t',header = T)
#allpair2=allpair2[allpair2$symbolnode1%in%c('MIAT'),]
feature=unique(c(allpair1$symbolnode1,allpair1$symbolnode2))

#phenotype
phenotype=read.table('~/repo/phenotypes_euro.txt',sep = ',',header = T,row.names = 1)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
phenotype=phenotype[phenotype$AtrialRhythm!='Yes AF/Sinus Rhythm',]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR')
rt=read.table('~/repo/rt.txt',sep = '\t',header = T,row.names = 1)
row.names(rt)=rt$name
#training set
train=rt[feature,row.names(phenotype)]

#get the matrix and sample information of GSE41177
exprSet=read.table("~/repo/exprs_gse41177.txt",sep = '\t',header = T)
sample=read.table("~/repo/sample_gse41177.txt",sep = '\t',header =T)
#training set
test=exprSet[feature,rownames(sample)]

#merge training set with test set
MERGE=cbind(train,test)
MERGE=na.omit(MERGE)

#remove batch effect 
library(sva)
combat_Data <- ComBat(dat = as.matrix(MERGE), batch = c(rep(1,ncol(train)),rep(2,ncol(test))))#not fit the interest in the model
MERGE=as.data.frame(t(MERGE))
MERGE$type=as.factor(c(phenotype$AtrialRhythm,sample$type))
colnames(MERGE)=str_replace(colnames(MERGE),'-','_')
TRAIN=MERGE[1:ncol(train),];
write.csv(TRAIN,file = '~/repo/TRAIN_LINC00964.csv',row.names = F)#use to condcuted 6 fold cross-validation of Random forest in Weka solfware
TEST=MERGE[(ncol(train)+1):nrow(MERGE),];


library(randomForest)
set.seed(54321)
rf_classifier <- randomForest(type ~ ., data=TRAIN, importance=TRUE,
                        proximity=TRUE)



# validation set assessment : ROC curves and AUC
prediction_for_roc_curve <- predict(rf_classifier,TEST[,-ncol(TEST)],type="prob")
library(pROC)
nb1.roc <- roc(response = as.numeric(TEST$type)-1,predictor = prediction_for_roc_curve[,1],levels = c(1, 0))
nb1.roc$auc#0.7894
plot.roc(nb1.roc,
         print.thres = TRUE,print.auc=TRUE, col="red",
         print.thres.best.method = "youden")



#using the same method to investigate the diagnostic role of MIAT's related ceRNA pairs