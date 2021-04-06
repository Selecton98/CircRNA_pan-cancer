load("IDCSCmatrix.Rda")
load("IDCSCmeta.Rda")
matrix3$circRNA<-NULL
#################################################
###loading output of WGCNA
node3q7<-read.table("./matrix3q7/nodes.txt",sep="\t",header=T)
node3q8<-read.table("./matrix3q8/nodes.txt",sep="\t",header=T)
node3q9<-read.table("./matrix3q9/nodes.txt",sep="\t",header=T)
node3q9orangered3<-node3q9[which(node3q9$nodeAttr.nodesPresent...=="orangered3"),]
node3q9firebrick4<-node3q9[which(node3q9$nodeAttr.nodesPresent...=="firebrick4"),]
node3q9salmon2<-node3q9[which(node3q9$nodeAttr.nodesPresent...=="salmon2"),]
node3q8coral4<-node3q8[which(node3q8$nodeAttr.nodesPresent...=="coral4"),]
node3q8sienna2<-node3q8[which(node3q8$nodeAttr.nodesPresent...=="sienna2"),]
node3q8darkolivegreen2<-node3q8[which(node3q8$nodeAttr.nodesPresent...=="darkolivegreen2"),]
node3q7whitesmoke<-node3q7[which(node3q7$nodeAttr.nodesPresent...=="whitesmoke"),]
node3q7deeppink2<-node3q7[which(node3q7$nodeAttr.nodesPresent...=="deeppink2"),]
orangered3<-colSums(matrix3[node3q9orangered3$nodeName,])
firebrick4<-colSums(matrix3[node3q9firebrick4$nodeName,])
salmon2<-colSums(matrix3[node3q9salmon2$nodeName,])
coral4<-colSums(matrix3[node3q8coral4$nodeName,])
sienna2<-colSums(matrix3[node3q8sienna2$nodeName,])
darkolivegreen2<-colSums(matrix3[node3q8darkolivegreen2$nodeName,])
whitesmoke<-colSums(matrix3[node3q7whitesmoke$nodeName,])
deeppink2<-colSums(matrix3[node3q7deeppink2$nodeName,])
metagene<-rbind(orangered3,firebrick4,salmon2,coral4,sienna2,darkolivegreen2,whitesmoke,deeppink2)
allcircs<-rbind(node3q9orangered3,node3q9firebrick4,node3q9salmon2,node3q8coral4,node3q8sienna2,node3q8darkolivegreen2,node3q7whitesmoke,node3q7deeppink2)
matrix3a<-matrix3[as.vector(allcircs$nodeName),]
library(pheatmap)
Samples<-as.data.frame(colnames(matrix3a))
colnames(Samples)<-"Samples"
meta5<-merge(Samples,meta4,all=F)
anno<-meta5
rownames(anno)<-anno$Samples
anno$Samples<-NULL
anno2<-allcircs[,c(1,3)]
rownames(anno2)<-anno2[,1]
anno2$nodeName<-NULL
ann_colors = list(nodeAttr.nodesPresent...=c(orangered3="orangered3",firebrick4="firebrick4",salmon2="salmon2",coral4="coral4",sienna2="sienna2",darkolivegreen2="darkolivegreen2",whitesmoke="lightgrey",deeppink2="deeppink2"),Type=c(cancer="tomato",normal="dodgerblue"))
matrix3a<-as.data.frame(t(matrix3a))
pheatmap(matrix3a,scale='column',annotation_row=anno,annotation_col=anno2,annotation_colors =ann_colors,cluster_col=T,cluster_row=T,color = colorRampPalette(c("navy","white", "red" ))(100),cellwidth=3,cellheight=1,fontsize =1) 
################################################
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