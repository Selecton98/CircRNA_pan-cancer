load("EXORBASEmatrix.Rda")
ALL[ALL< 2] <- 0
rownames(ALL)<-ALL$circ
ALL$circ<-NULL
ALL<-ALL[which(rowSums(ALL)>0),]
ALL<-ALL[,which(colSums(ALL)>0)]
ALL2<-ALL
ALL2$circ<-rownames(ALL2)
#######################################################
###figure 3.1
cleancan<-as.data.frame(table(rowSums(as.matrix(ALL[1:21])>0)))
cleannor<-as.data.frame(table(rowSums(as.matrix(ALL[22:53])>0)))
cleancan<-cleancan[-1,]
cleannor<-cleannor[-1,]
cleancan$Var1<-as.character(cleancan$Var1)
cleannor$Var1<-as.character(cleannor$Var1)
cleancan$Var1<-as.numeric(cleancan$Var1)
cleannor$Var1<-as.numeric(cleannor$Var1)
cleancan$tag<-"cancer"
cleannor$tag<-"normal"
cleancan$Var1p<-cleancan$Var1/21
cleannor$Var1p<-cleannor$Var1/32
cleancannor<-rbind(cleancan,cleannor)
cleancannor$logfreq<-log10(cleancannor$Freq)
ggplot(data=cleancannor, mapping=aes(x=Var1p, y=logfreq,group=tag)) +
  geom_line(aes(color=tag))+
  scale_x_continuous(name="Proportion of samples haboring the circRNA",breaks=c(0,0.25,0.50,0.75,1.00))+
  scale_y_continuous(name="log10(Types of the circRNA)",breaks=c(0,2,2.5,3,3.5,4))+
  theme_classic() 
#######################################################
###figure3.2
samplesumcan<- as.data.frame(colSums(ALL[1:21]))
samplesumnor<- as.data.frame(colSums(ALL[22:53]))
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
ggplot(data=samplesumall, aes(x=tag, y=Var1,fill= tag))+
  geom_violin(alpha=0.8,width=1)+
  xlab("tag")+
  ylab("Total count of circRNA in each sample")+
  theme_classic()
#######################################################
###figure 3.3
amplesumcan<- as.data.frame(rowSums(ALL[1:21]))
amplesumnor<- as.data.frame(rowSums(ALL[22:53]))
amplesumcan$tag<-"cancer"
amplesumnor$tag<-"normal"
colnames(amplesumcan)<-c("Var1","tag")
colnames(amplesumnor)<-c("Var1","tag")
amplesumcan<-amplesumcan[which(amplesumcan$Var1>0),]
amplesumnor<-amplesumnor[which(amplesumnor$Var1>0),]
amplesumall<-rbind(amplesumcan,amplesumnor)
amplesumall$Var1<-as.character(amplesumall$Var1)
amplesumall$Var1<-as.numeric(amplesumall$Var1)
amplesumall$Var1p<-log10(amplesumall$Var1)
ggplot(data=amplesumall, aes(x=tag, y=Var1p,fill= tag))+
  geom_violin(alpha=0.8,width=1)+  
  xlab("tag")+
  ylab("log10(Total count of a specific circRNA in all samples)")+
  theme_classic()
##############################################
###wilcox test
colnames(ALL)<-gsub("clean","",colnames(ALL))
colnames(ALL)<-gsub("CIRCexplorer_circ.txt","",colnames(ALL))
ALLt<-as.data.frame(t(ALL))
rownames(ALLt)
Y<-as.data.frame(c(rep("HCC",21),rep("NP",32)))
colnames(Y)<-"Y"
ALLt<-cbind(Y,ALLt)
mixed<-read.csv(file="lasso-mixed-gene.csv")
nor<-read.csv(file="lasso-matrix3q93-gene.csv")
together<-rbind(mixed,nor)
together2<-intersect(together$circRNA,colnames(ALLt))
together2<-c("Y",together2)
ALLfor<-ALLt[,as.vector(together2)]
library(ggpubr)
library(RColorBrewer)
my_comparisons <- list(c("HCC","NP"))
for(i in 2:ncol(ALLfor)){
  df1<-ALLfor[,c(1,i)]
  name<-colnames(df1)[2]
  p<-ggviolin(df1, x="Y", y=name, fill = "Y", 
              palette = c("tomato","dodgerblue"))+
    stat_compare_means(comparisons = my_comparisons)+#label这里表示选择显著性标记（星号） 
    stat_compare_means(label.y = 50)
  p
  myfilename <- paste(name,"_violin.pdf",sep="")
  ggsave(filename = myfilename,width = 7, height = 7)
}
unique<-read.csv(file="cancerspecific-cancercirc5.csv")
unique2<-intersect(unique$circRNA,colnames(ALLt))
unique2<-c("Y",unique2)
ALLfor<-ALLt[,as.vector(unique2)]
library(ggpubr)
library(RColorBrewer)
my_comparisons <- list(c("HCC","NP"))
for(i in 2:ncol(ALLfor)){
  df1<-ALLfor[,c(1,i)]
  name<-colnames(df1)[2]
  p<-ggviolin(df1, x="Y", y=name, fill = "Y", 
              palette = c("tomato","dodgerblue"))+
    stat_compare_means(comparisons = my_comparisons)+#label这里表示选择显著性标记（星号） 
    stat_compare_means(label.y = 50)
  p
  myfilename <- paste(name,"_violin.pdf",sep="")
  ggsave(filename = myfilename,width = 7, height = 7)
}
