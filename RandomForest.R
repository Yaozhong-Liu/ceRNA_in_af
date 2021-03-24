#Trans=read.csv("G:/生信分析/Integrative omics data/BMC 10.13/REM.csv")
Trans=read.csv("G:/生信分析/Integrative omics data/GeneMeta/REM.csv")
allpair1=read.table('G:/生信分析/ceRNA/2020.12/allpair_turquoise.txt',sep = '\t',header = T)
allpair1=allpair1[allpair1$symbolnode1%in%c('LINC00964'),]
lyz1=Trans[Trans$X%in%c(allpair1$symbolnode1,allpair1$symbolnode2),]

allpair2=read.table('G:/生信分析/ceRNA/2020.12/allpair_tan.txt',sep = '\t',header = T)
allpair2=allpair2[allpair2$symbolnode1%in%c('MIAT'),]
#lyz2=Trans[Trans$X%in%c(allpair2$symbolnode1,allpair2$symbolnode2),]
#mean(lyz2$zSco_Ex_8)


feature=unique(c(allpair1$symbolnode1,allpair1$symbolnode2))

#####-------------------------------------------------------------------------------------------------------------------------------------------
phenotype=read.table('G:/生信分析/ceRNA/2020.12/GSE69890_RAW/phenotypes_euro.txt',sep = ',',header = T,row.names = 1)
phenotype=as.data.frame(t(phenotype))
phenotype=phenotype[,c(1,30)]
phenotype=phenotype[phenotype$AtrialRhythm!='Yes AF/Sinus Rhythm',]
library(stringr)
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/AF Rhythm','AF')
phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'No AF/Sinus Rhythm','SR')
#phenotype$AtrialRhythm=str_replace(phenotype$AtrialRhythm,'Yes AF/Sinus Rhythm','AF')

rt=read.table('G:/生信分析/ceRNA/2020.12/rt.txt',sep = '\t',header = T,row.names = 1)
row.names(rt)=rt$name
train=rt[feature,row.names(phenotype)]




#----------
#exprSet=read.csv("G:/生信分析/Meta+machine learning/merged/exprs.csv",header = T,row.names = 1)
#sample=read.csv("G:/生信分析/Meta+machine learning/merged/sample.csv",header =T,row.names = 1)
#sample=sample[1:130,]
#sample=sample[sample$GSE=='GSE14975',]
exprSet=read.csv("G:/生信分析/AF dataset/AF+mRNA/GSE41177/exprs.csv",header = T,row.names = 1)
sample=read.csv("G:/生信分析/AF dataset/AF+mRNA/GSE41177/sample.csv",header =T,row.names = 1)

table(sample$type)
#sample=sample[sample$origin!=8,]
#gse2240=read.table("G:/生信分析/AF dataset/AF+mRNA/GSE2240_RAW/exprs.txt",sep = '\t',header = T,row.names = 1)
#lyz=gse2240[row.names(gse2240)%in%c(allpair1$symbolnode1,allpair1$symbolnode2),]
#library(limma)
#class <- sample[sample$origin==8,]$type
#design <- model.matrix(~-1+factor(class))
#colnames(design) <- c("AF","SR")
#contrast.matrix=makeContrasts(contrasts = "AF-SR",levels = design)
#fit <- lmFit(gse2240,design)
#fit1=contrasts.fit(fit,contrast.matrix)
#fit2 <- eBayes(fit1)
#allDiff=topTable(fit2,adjust='fdr',coef="AF-SR",number=200000)
#allDiff[row.names(allDiff)%in%c(allpair2$symbolnode1,allpair2$symbolnode2),]
#WEKA=exprSet[unique(c(allpair1$symbolnode1,allpair1$symbolnode2,allpair2$symbolnode1,allpair2$symbolnode2)),]
#WEKA["type",]=sample[1:130,]$type
#WEKA=as.data.frame(t(WEKA))
#write.csv(WEKA,file = 'G:/生信分析/ceRNA/2020.12/WEKA.csv',row.names = F)
test=exprSet[feature,rownames(sample)]

#test=exprSet[unique(c(allpair2$symbolnode1,allpair2$symbolnode2)),rownames(sample)]
MERGE=cbind(train,test)
MERGE=na.omit(MERGE)

#train=exprSet[unique(c(allpair1$symbolnode1,allpair1$symbolnode2,allpair2$symbolnode1,allpair2$symbolnode2)),]
#gse2240=read.table("G:/生信分析/AF dataset/AF+mRNA/GSE2240_RAW/exprs.txt",sep = '\t',header = T,row.names = 1)
#test=gse2240[row.names(train),]
#MERGE=cbind(train,test)
#MERGE=na.omit(MERGE)



#a=intersect(rownames(rt),row.names(exprSet))
#MERGE=cbind(rt[a,row.names(phenotype)],exprSet[a,row.names(sample)])
#MERGE=MERGE[unique(c(allpair1$symbolnode1,allpair1$symbolnode2,allpair2$symbolnode1,allpair2$symbolnode2)),]

#去除批次效应
#Combat去除batch效应
library(sva)
#modcombat = model.matrix(~as.factor(sample$type)) 
combat_Data <- ComBat(dat = as.matrix(MERGE), batch = c(rep(8,ncol(train)),sample$GSE))#,mod = factor(c(phenotype$AtrialRhythm,sample$type))
MERGE=as.data.frame(combat_Data)#,
#datExpr=t(MERGE)
#D=dist(datExpr,diag = FALSE,upper = FALSE)
#hc=hclust(D)
#plot(hc)
#rect.hclust(hc,2)
MERGE=as.data.frame(t(MERGE))
MERGE$type=c(phenotype$AtrialRhythm,sample$type)
MERGE$type=as.factor(MERGE$type)
colnames(MERGE)=str_replace(colnames(MERGE),'-','_')
TRAIN=MERGE[1:ncol(train),];#write.csv(TRAIN,file = 'G:/生信分析/ceRNA/2020.12/TRAIN_LINC00964.csv',row.names = F)
TEST=MERGE[(ncol(train)+1):nrow(MERGE),];#write.csv(TEST,file = 'G:/生信分析/ceRNA/2020.12/TEST_LINC00964.csv',row.names = F)
#MERGE$type=sample$type
#MERGE$type=as.factor(MERGE$type)
#TRAIN=MERGE[1:130,];#
#TEST=MERGE[131:160,];#


library(randomForest)
##data(iris)
set.seed(54321)
rf_classifier <- randomForest(type ~ ., data=TRAIN, importance=TRUE,
                        proximity=TRUE)

print(rf_classifier)
varImpPlot(rf_classifier)


# Validation set assessment #1: looking at confusion matrix
prediction_for_table <- predict(rf_classifier,TEST[,-ncol(TEST)])
table(observed=TEST[,ncol(TEST)],predicted=prediction_for_table)

# Validation set assessment #2: ROC curves and AUC
# Calculate the probability of new observations belonging to each class
# prediction_for_roc_curve will be a matrix with dimensions data_set_size x number_of_classes
prediction_for_roc_curve <- predict(rf_classifier,TEST[,-ncol(TEST)],type="prob")
library(pROC)
nb1.roc <- roc(response = as.numeric(TEST$type)-1,predictor = prediction_for_roc_curve[,1],levels = c(1, 0))
nb1.roc$auc#0.7894
plot.roc(nb1.roc,
         print.thres = TRUE,print.auc=TRUE, col="red",
         print.thres.best.method = "youden")
pred1 <- ifelse(test = prediction_for_roc_curve[,1] >= 0.552, 
                    yes = "AF", #... assign the positive class (No)
                    no = "SR") #... assign the negative class (Yes)
pred1 <- as.factor(pred1)
table(observed=TEST[,ncol(TEST)],predicted=pred1)

?plot.roc

###NB
#library(e1071)
#nb1 <- naiveBayes(type ~ ., data = TRAIN)
#print(nb1)
#nb1.pred <- predict(nb1, newdata = TEST, type = 'class')
#table(true = TEST$type, predicted = nb1.pred)
#nb1.pred.prob <- predict(nb1, newdata = TEST, type = "raw")
#nb1.pred.prob
#library(pROC)
#nb1.roc <- roc(response = as.numeric(TEST$type),predictor = nb1.pred.prob[,1],levels = c(2, 1))
#nb1.roc$auc
#plot.roc(nb1.roc,
#         print.thres = TRUE,
#         print.thres.best.method = "youden")






#diffExpLevel=t(TEST[,-174])
#diffExpLevel <- exprSet[unique(c(allpair1$symbolnode1,allpair1$symbolnode2,allpair2$symbolnode1,allpair2$symbolnode2)),row.names(sample[order(sample$type),])]
diffExpLevel=na.omit(diffExpLevel)
library(pheatmap)
annotation <- sample[,c(1,7)]
#绘制热图=====================================
pheatmap(diffExpLevel, #表达数据
         cluster_rows = T,#行聚类
         cluster_cols = T,#列聚类
         annotation_col =annotation, #样本分类数据
         annotation_legend=TRUE, # 显示样本分类
         show_rownames = F,# 显示行名
         show_colnames = T,# 显示列名
         scale = "row", #对行标准化
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100) # 热图基准颜色
)
#dev.off()


