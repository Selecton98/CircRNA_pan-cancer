load("IDCSCmatrix.Rda")
load("IDCSCmeta.Rda")
rownames(matrix3)<-matrix3$circRNA
matrix3$circRNA<-NULL
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.6)
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.7)
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.8)
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.9)
matrix3q7<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=4&rowSums(as.matrix(matrix3)> 0)<7),]
matrix3q8<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=7&rowSums(as.matrix(matrix3)> 0)<20),]
matrix3q9<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=20),]
################################################################
###1 ML data preparation
matrix3t<-as.data.frame(t(matrix3q9)) #according to which quantile to be investigated
Samples<-as.data.frame(rownames(matrix3t))
colnames(Samples)<-"Samples"
meta5<-merge(Samples,meta4,all=F)
tag<-gsub("cancer","4",meta5$Type)
tag<-gsub("normal","2",tag)
matrix3t<-cbind(tag,matrix3t)
# partitioning into training (70%) and validation (30%)
set.seed(42)
# randomly sample 70% of the row IDs for training; the remaining 30% serve as validation
train.rows <- sample(rownames(matrix3t), dim(matrix3t)[1]*0.7)
# collect all the columns with training row ID into training set:
matrix3t.train <- matrix3t[train.rows, ]
# assign row IDs that are not already in the training set, into validation
valid.rows <- setdiff(rownames(matrix3t), train.rows)
matrix3t.valid <- matrix3t[valid.rows, ]
##################################################################
###2 feature selection by lasso
set.seed(42) 
library(glmnet)
library(pROC)
lasso.results <- c()
for (j in 1:50) {
  s <- matrix3t.train[sample(1:nrow(matrix3t.train),ceiling(dim(matrix3t.train)[1]*0.5),replace=F),]
  ts <- apply(s,2,as.numeric)
  y <- as.matrix(ts[,1])
  x <- as.matrix(ts[,-1])
  cv.fit <- cv.glmnet(x,y,family="binomial", type.measure = "auc", nfolds=5)
  co<-coef(cv.fit,s="lambda.1se")
  name <- rownames(co)[co[,1]!=0]
  lasso.results <- c(lasso.results, name)
  pred = predict(cv.fit, newx = x, type = 'response',s ="lambda.min")
  roc<-roc(y,pred)
  pdf(paste(j,".ROC01.pdf",sep=""))
  plot.roc(roc,col="red",print.auc =TRUE,print.auc.col = "darkgreen",auc.polygon = TRUE,auc.polygon.col = "pink")
 dev.off() }
plot(cv.fit)
f1 = glmnet(x, y, family="binomial",type.measure = "auc", nfolds=5)
plot(f1, xvar="lambda", label=TRUE)
freq.lasso.results <- as.data.frame(table(lasso.results))
freq.lasso.results<-freq.lasso.results[-1,]
write.csv(freq.lasso.results,"matrix3q7lasso.csv")
####################################################################
###3 feature selection by random forest
library(varSelRF)
facy <- factor(matrix3t.train[,1])
x <- matrix3t.train[,-1]
step=varSelRF(x, facy, c.sd = 1, mtryFactor = 1, ntree = 500,
              ntreeIterat = 500, vars.drop.num = NULL, vars.drop.frac = 0.1,
              whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE,
              returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE)
select.history <- step$selec.history
selected.vars<-step$selected.vars
selected.vars2<-as.data.frame(selected.vars)
write.table(select.history,"20210201select.history.txt", row.names = F, quote = F)
write.csv(selected.vars2,"2020201RF.csv")
######################################################################
###4 combine lasso and random forest
filter.lasso<-freq.lasso.results[which(freq.lasso.results$Freq>15),]
filter.lasso<-filter.lasso[,-2]
lasso<-as.vector(filter.lasso)
RF<-selected.vars
library(VennDiagram)
input<-list(RF, lasso)
Table<-calculate.overlap(input)
venn<-venn.diagram(input,NULL, main.cex = 3,
                   category = c("Random forest","Lasso"),fill = c( "orange","lightgrey"),
                   cat.col= c("orange","lightgrey"),   imagetype = "tiff",  main.fontfamily="serif") #, filename = "Vennup.tif"
grid.draw(venn)
sect<-as.vector(filter.lasso)
sect[length(sect)+1]<-"tag"
matrix3t.train3<-matrix3t.train[,sect]
matrix3t.valid3<-matrix3t.valid[,sect]
matrix3t3<-matrix3t[,sect]
colnames(matrix3t.train3)<-gsub("tag","Y",colnames(matrix3t.train3))
colnames(matrix3t.valid3)<-gsub("tag","Y",colnames(matrix3t.valid3))
colnames(matrix3t3)<-gsub("tag","Y",colnames(matrix3t3))
###################################################################
###5 SVM diganostic model construction
set.seed(42)
library(pROC)
library(caret)
library(e1071)
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)
matrix3t.train3$Y<-gsub("4","cancer",matrix3t.train3$Y)   
matrix3t.train3$Y<-gsub("2","normal",matrix3t.train3$Y) 
matrix3t.valid3$Y<-gsub("4","cancer",matrix3t.valid3$Y)   
matrix3t.valid3$Y<-gsub("2","normal",matrix3t.valid3$Y) 
svmFit <- train(Y ~ ., data=matrix3t.train3, 
                method = "svmRadial", 
                trControl = fitControl, 
                #          preProc = c("center", "scale"),
                tuneLength = 10,
                metric = "ROC")
train.pred <- predict(svmFit,matrix3t.train3)
validation.pred <- predict(svmFit,matrix3t.valid3)
tpred <- predict(svmFit,matrix3t.train3, type = 'prob')
vpred <- predict(svmFit,matrix3t.valid3, type = 'prob')
train.con <- confusionMatrix(train.pred, factor(matrix3t.train3$Y))
valid.con <- confusionMatrix(validation.pred, factor(matrix3t.valid3$Y))
#testing
testing.pred <- predict(svmFit,HCCNP3)
testing.con <- confusionMatrix(testing.pred, factor(HCCNP3$Y))
spred <- predict(svmFit,HCCNP3, type = 'prob')
train.con
valid.con 
testing.con
##################################################################
###6 plot ROC for SVM
z<-gsub("cancer","4",matrix3t.train3$Y) 
z<-gsub("normal","2",z)   
q<-gsub("cancer","4",matrix3t.valid3$Y) 
q<-gsub("normal","2",q) 
w<-gsub("cancer","4",HCCNP3$Y) 
w<-gsub("normal","2",w) 
roc1<-roc(q,vpred[,1])
plot.roc(roc1,col="red",print.auc =TRUE,print.auc.col = "red")
roc2<-roc(z,tpred[,1])
plot.roc(roc2,col="dodgerblue",print.auc =TRUE,print.auc.col = "dodgerblue",add=TRUE)
roc3<-roc(w,spred[,1])
plot.roc(roc3,col="darkgreen",print.auc =TRUE,print.auc.col = "darkgreen",add=TRUE)
library(pheatmap)
matrixtsne<-matrix3t3
tg<-gsub("2","normal",matrixtsne$Y) 
tg<-gsub("4","cancer",tg) 
anno<-as.data.frame(rownames(matrixtsne),tg)
anno$Label<-rownames(anno)
rownames(anno)<-anno$`rownames(matrixtsne)`
anno$`rownames(matrixtsne)`<-NULL
matrixtsne$Y<-NULL
ann_colors = list(Label=c(cancer="tomato",normal="dodgerblue"))
pheatmap(matrixtsne,scale='column',annotation_row=anno,annotation_colors =ann_colors,cluster_col=T,cluster_row=T,color = colorRampPalette(c("navy","white", "red"))(100),cellwidth=10,cellheight=1,fontsize =1) 