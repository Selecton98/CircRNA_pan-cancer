load("IDCSCmatrix.Rda")
load("IDCSCmeta.Rda")
rownames(matrix3)<-matrix3$circRNA
matrix3$circRNA<-NULL
load("circgeneconvert.Rda")
matrix3$circRNA<-rownames(matrix3)
matrix3anno<-merge(circgene3,matrix3,all=F) 
meta5<-meta4[,c(1,2)]
meta6<-meta5[which(meta5$Type=="cancer"),]
meta7<-meta5[which(meta5$Type=="normal"),]
meta6<-as.vector(meta6$Samples)
meta7<-as.vector(meta7$Samples)
meta8<-c("circRNA",meta6)
meta9<-c("circRNA",meta7)
meta6<-c("gene_name",meta6)
meta7<-c("gene_name",meta7)
matrix3canceranno<-matrix3anno[,meta6]
matrix3normalanno<-matrix3anno[,meta7]
matrix3cancercirc<-matrix3[,meta8]
matrix3normalcirc<-matrix3[,meta9]
matrix3cancercirc<-matrix3cancercirc[which(rowSums(matrix3cancercirc[2:266])>0),]
matrix3normalcirc<-matrix3normalcirc[which(rowSums(matrix3normalcirc[2:320])>0),]
cleancan<-as.data.frame(table(rowSums(as.matrix(matrix3cancercirc[2:266])>0)))
cleannor<-as.data.frame(table(rowSums(as.matrix(matrix3normalcirc[2:320])>0)))
samplesumcan<- as.data.frame(colSums(matrix3cancercirc[2:266]))
samplesumnor<- as.data.frame(colSums(matrix3normalcirc[2:320]))
amplesumcan<- as.data.frame(rowSums(matrix3cancercirc[2:266]))
amplesumnor<- as.data.frame(rowSums(matrix3normalcirc[2:320]))
cleancan$Var1<-as.character(cleancan$Var1)
cleannor$Var1<-as.character(cleannor$Var1)
cleancan$Var1<-as.numeric(cleancan$Var1)
cleannor$Var1<-as.numeric(cleannor$Var1)
cleancan$tag<-"cancer"
cleannor$tag<-"normal"
cleancan$Var1p<-cleancan$Var1/265
cleannor$Var1p<-cleannor$Var1/319
cleancannor<-rbind(cleancan,cleannor)
cleancannor$logfreq<-log10(cleancannor$Freq)
##############################################
###figure 1.0
meta4can<-meta4[which(meta4$Type=="cancer"),]
sampleident<-as.data.frame(table(meta4can$Site))
sampleident<-sampleident[which(sampleident$Freq!=0),]
sum(sampleident$Freq)
sampleident$Var1<-paste(sampleident$Var1,"neoplasms",sep="_")
meta4nor<-meta4[which(meta4$Type=="normal"),]
sampleidentn<-as.data.frame(table(meta4nor$Site))
sampleidentn<-sampleidentn[which(sampleidentn$Freq!=0),]
sum(sampleidentn$Freq)
sampleident<-rbind(sampleident,sampleidentn)
ggplot(sampleident,aes(x= reorder(Var1,Freq), y=Freq,fill=Var1)) +                                                                                             geom_bar(stat = "identity") +
xlab("Sites of normal tissues") +
coord_flip()+
scale_y_continuous(name="Numbers of samples",expand=c(0,0))+
theme_classic()
###############################################
###figure1.1
summary(cleancan$Var1p)
summary(cleannor$Var1p)
ggplot(data=cleancannor, mapping=aes(x=Var1p, y=logfreq,group=tag)) +
geom_line(aes(color=tag))+
scale_x_continuous(name="Proportion of samples haboring the circRNA",breaks=c(0,0.25,0.50,0.75,0.82,0.98,1.00))+
scale_y_continuous(name="log10(Types of the circRNA)")+
theme_classic() 
###figure1.2
samplesumcan$tag<-"cancer"
samplesumnor$tag<-"normal"
colnames(samplesumcan)<-c("Var1","tag")
colnames(samplesumnor)<-c("Var1","tag")
samplesumall<-rbind(samplesumcan,samplesumnor)
samplesumall$Var1<-as.character(samplesumall$Var1)
samplesumall$Var1<-as.numeric(samplesumall$Var1)
ggplot(data=samplesumall,aes(x=Var1,fill=factor(tag),alpha=0.4))+
geom_density(stat="density")+
scale_x_continuous(name="Total count of circRNA in individual samples",expand=c(0,0))+
scale_y_continuous(name="Density",expand=c(0,0))+
theme_classic()
##################################################
###figure1.3
ggplot(data=samplesumall, aes(x=tag, y=Var1,fill= tag))+
geom_violin(alpha=0.8,width=1)+
xlab("tag")+
ylab("Total count of circRNA in each sample")+
theme_classic()
####################################################
###figure1.4
amplesumcan$tag<-"cancer"
amplesumnor$tag<-"normal"
colnames(amplesumcan)<-c("Var1","tag")
colnames(amplesumnor)<-c("Var1","tag")
amplesumall<-rbind(amplesumcan,amplesumnor)
amplesumall$Var1<-as.character(amplesumall$Var1)
amplesumall$Var1<-as.numeric(amplesumall$Var1)
summary(amplesumall$Var1)
amplesumall$Var1p<-log10(amplesumall$Var1)
ggplot(data=amplesumall, aes(x=tag, y=Var1p,fill= tag))+
geom_violin(alpha=0.8,width=1)+  
xlab("tag")+
ylab("log10(Total count of a specific circRNA in all samples)")+
theme_classic()
matrix3canceranno<-matrix3canceranno[which(rowSums(matrix3canceranno[2:266])>0),]
matrix3normalanno<-matrix3normalanno[which(rowSums(matrix3normalanno[2:320])>0),]
matrix3canceru<-matrix3canceranno %>% group_by(gene_name) %>% summarise_each(funs(sum)) 
matrix3normalu<-matrix3normalanno %>% group_by(gene_name) %>% summarise_each(funs(sum)) 
tablecircgenecancer<-as.data.frame(table(matrix3canceranno$gene_name))
tablecircgenenormal<-as.data.frame(table(matrix3normalanno$gene_name))
########################################################
### figure2.1
library(ggplot2)
tablecancerp<-tablecircgenecancer[which(tablecircgenecancer$Freq>60),]
howmany<-unique(matrix3cancer2$circRNA)
tablecancerp$Freqp<-tablecancerp$Freq/64180
ggplot(tablecancerp,aes(x= reorder(Var1,Freq), y=Freq,fill=Freq)) +                                                                                            
geom_bar(stat = "identity") +
coord_flip()+
scale_fill_gradient(low = "pink", high = "red")+
xlab("Highly back-spliced host genes in cancer") +
scale_y_continuous(name="Types of circRNA",expand=c(0,0))+
theme_classic()
ggplot(tablecancerp,aes(x= reorder(Var1,Freqp), y=Freqp,fill=Freqp)) +                                                                                      geom_bar(stat = "identity") +
coord_flip()+
scale_fill_gradient(low = "pink", high = "red")+
xlab("Highly back-spliced host genes in cancer") +
scale_y_continuous(name="Types of circRNA",expand=c(0,0))+
theme_classic()
###################################################
###figure 2.7
viocVar1<-as.data.frame(tablecircgenecancer$Freq)
viocVar1$tag<-"cancer"
colnames(viocVar1)<-c("Var1","tag")
vionVar1<-as.data.frame(tablecircgenenormal$Freq)
vionVar1$tag<-"normal"
colnames(vionVar1)<-c("Var1","tag")
vio<-rbind(vionVar1,viocVar1)
vio$Var1p<-log10(vio$Var1)
ggplot(data=vio, aes(x=tag, y=Var1,fill= tag))+
xlab("tag")+
ylab("Total count of circRNA sharing a specific host gene")+
geom_violin(alpha=0.8,width=1)+  
theme_classic()
ggplot(data=vio, aes(x=tag, y=Var1p,fill= tag))+
xlab("tag")+
ylab("log10(Total count of circRNA sharing a specific host gene)")+
geom_violin(alpha=0.8,width=1)+  
theme_classic()
######################################################
###figure2.2
tablenormalp<-tablecircgenenormal[which(tablecircgenenormal$Freq>180),]
howmany2<-unique(matrix3normal2$circRNA)
tablenormalp$Freqp<-tablenormalp$Freq/200821
ggplot(tablenormalp,aes(x= reorder(Var1,Freq), y=Freq,fill=Freq)) +                                                                                           
geom_bar(stat = "identity") +
coord_flip()+
scale_fill_gradient(low = "dodgerblue", high = "navy")+
xlab("Highly back-spliced host genes in normal") +
scale_y_continuous(name="Types of circRNA",expand=c(0,0))+
theme_classic()
ggplot(tablenormalp,aes(x= reorder(Var1,Freqp), y=Freqp,fill=Freqp)) +                                                                                     geom_bar(stat = "identity") +
coord_flip()+
scale_fill_gradient(low = "dodgerblue", high = "navy")+
xlab("Highly back-spliced host genes in normal") +
scale_y_continuous(name="Types of circRNA",expand=c(0,0))+
theme_classic()
write.csv(tablecircgenenormal,"tablecircgenenormal.csv")
write.csv(tablecircgenecancer,"tablecircgenecancer.csv")
write.csv(tablenormalp,"tablenormalp.csv")
write.csv(tablecancerp,"tablecancerp.csv")
######################################################
###figure2.3
library(VennDiagram)
input<-list(tablecircgenenormal$Var1,tablecircgenecancer$Var1)
Table<-calculate.overlap(input)
venn<-venn.diagram(input,NULL, main.cex = 3,
category = c("Host_gene_normal","Host_gene_cancer"),fill = c( "dodgerblue","tomato"),
cat.col= c("dodgerblue","tomato"),   imagetype = "tiff",  main.fontfamily="serif") #, filename = "Vennup.tif"
grid.draw(venn)
commongene<-intersect(tablecircgenenormal$Var1,tablecircgenecancer$Var1)
cancergene<-setdiff(tablecircgenecancer$Var1,commongene)
normalgene<-setdiff(tablecircgenenormal$Var1,commongene)
cancergene<-as.data.frame(cancergene)
normalgene<-as.data.frame(normalgene)
write.csv(cancergene,"cancergene.csv")
write.csv(normalgene,"normalgene.csv")
#######################################################
###figure2.8
save<-vector()
for (j in 1:100) {
i<-j/100
tablecancerq<-tablecircgenecancer[which(tablecircgenecancer$Freq>=quantile(tablecircgenecancer$Freq,i)),]
tablenormalq<-tablecircgenenormal[which(tablecircgenenormal$Freq>=quantile(tablecircgenenormal$Freq,i)),]
commongene<-intersect(tablecancerq$Var1,tablenormalq$Var1)
save[j]<-length(commongene)}
changeoverlap<-data.frame(c(1:100),save)
colnames(changeoverlap)<-c("sequencer","overlap")
changeoverlap$quantilel<-changeoverlap$sequencer/100
ggplot(data=changeoverlap, mapping=aes(x=quantilel, y=overlap,color=overlap)) +
geom_line()+
scale_color_gradient()+
scale_x_continuous(name="Quantile of types of circRNAs",breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))+
scale_y_continuous(name="Number of overlapping host genes",breaks=c(0,50,100,200,300,400,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500))+theme_classic() 
########################################################
###figure1.5
library(VennDiagram)
input<-list(matrix3normalcirc$circRNA,matrix3cancercirc$circRNA)
Table<-calculate.overlap(input)
venn<-venn.diagram(input,NULL, main.cex = 3,
category = c("CircRNA_normal","CircRNA_cancer"),fill = c( "dodgerblue","tomato"),
cat.col= c("dodgerblue","tomato"),   imagetype = "tiff",  main.fontfamily="serif") #, filename = "Vennup.tif"
grid.draw(venn)
commoncirc<-intersect(matrix3normalcirc$circRNA,matrix3cancercirc$circRNA)
cancercirc<-setdiff(matrix3cancercirc$circRNA,commoncirc)
normalcirc<-setdiff(matrix3normalcirc$circRNA,commoncirc)
cancercirc<-as.data.frame(cancercirc)
normalcirc<-as.data.frame(normalcirc)
write.csv(cancercirc,"cancercirc.csv")
write.csv(normalcirc,"normalcirc.csv")
###########################################################
###figure2.4
library(ggplot2)
ggplot(data=cleancanunor, mapping=aes(x=Var1, y=logfreq,group=tag)) +
geom_line(aes(color=tag))+
theme_classic() +
scale_x_continuous(name="Number of samples haboring the host gene-circRNA")+
scale_y_continuous(name="log10(Types of host gene-circRNA)")
summary(cleancanu$Var1p)
summary(cleannoru$Var1p)
ggplot(data=cleancanunor, mapping=aes(x=Var1p, y=logfreq,group=tag)) +
geom_line(aes(color=tag))+
scale_x_continuous(name="Proportion of samples haboring the circRNA",breaks=c(0,0.25,0.50,0.75,0.85,0.98,1.00))+
scale_y_continuous(name="log10(Types of the circRNA)")+
theme_classic() 