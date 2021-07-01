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
###################################################################
###3 SVM diganostic model construction
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
###4 plot ROC for SVM
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
#################################################################
###loading output of LASSO
lasso9<-read.csv("matrix3q9lasso.csv",header=T)
lasso8<-read.csv("matrix3q8lasso.csv",header=T)
lasso7<-read.csv("matrix3q7lasso.csv",header=T)
lasso9$X<-NULL
lasso8$X<-NULL
lasso7$X<-NULL
lasso9<-lasso9[which(lasso9$Freq>10),]
lasso8<-lasso8[which(lasso8$Freq>15),]
lasso7<-lasso7[which(lasso7$Freq>15),]
library("ggplot2")
ggplot(lasso7,aes(x= reorder(lasso.results,Freq), y=Freq,fill=Freq)) +   
geom_bar(stat = "identity") +
coord_flip()+
scale_fill_gradient(low = "pink", high = "red")+
xlab("circRNAs") +
scale_y_continuous(name="Frequency",expand=c(0,0))+
theme_classic()
lasso9$Block<-"Top10"
lasso8$Block<-"Top1020"
lasso7$Block<-"Top2030"
lasso<-rbind(lasso9,lasso8,lasso7)
matrix3b<-matrix3[as.vector(lasso$lasso.results),]
library(pheatmap)
library(RColorBrewer)
anno2<-lasso[,c(1,3)]
rownames(anno2)<-anno2[,1]
anno2$lasso.results<-NULL
ann_colors = list(Block=c(Top10=brewer.pal(3,"Set3")[1],Top1020=brewer.pal(3,"Set3")[2],Top2030=brewer.pal(3,"Set3")[3]),Type=c(cancer="tomato",normal="dodgerblue"))
matrix3b<-as.data.frame(t(matrix3b))
pheatmap(matrix3b,scale='column',annotation_row=anno,annotation_col=anno2,annotation_colors =ann_colors,cluster_col=T,cluster_row=T,color = colorRampPalette(c("navy","white", "red" ))(100),cellwidth=3,cellheight=1,fontsize =1) 
write.csv(colnames(matrix3b),"lassotopcircRNA.csv")
