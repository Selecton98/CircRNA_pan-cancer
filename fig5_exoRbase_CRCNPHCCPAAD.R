library(dplyr)
library(tidyr)
library(DESeq2)
load(file="CRCNPHCCPAAD.Rda")
#####################################
rownames(big)<-big$circ
big$circ<-NULL
know2<-as.data.frame(t(big))

CRCname<-data.frame(colnames(CRC),"CRC")
HCCname<-data.frame(colnames(HCC),"HCC")
NPname<-data.frame(colnames(NP),"NP")
PAADname<-data.frame(colnames(PAAD),"PAAD")

colnames(CRCname)<-c("sample","Type")
CRCname<-CRCname[-1,]

colnames(PAADname)<-c("sample","Type")
PAADname<-PAADname[-1,]

colnames(NPname)<-c("sample","Type")
NPname<-NPname[-1,]

colnames(HCCname)<-c("sample","Type")
HCCname<-HCCname[-1,]

type<-rbind(CRCname,PAADname,NPname,HCCname)

rownames(type)<-type$sample
type$sample<-NULL

type<-type[rownames(know2),]
type<-as.data.frame(type)
colnames(type)<-"Type"

#############################################
##PCA analysis

library(ggplot2)
library(ggfortify)
library(devtools)
know2<-cbind(know2,type)
autoplot(prcomp(know2[,1:73304]), data=know2,label=FALSE,label.size=4,colour="Type")+theme_bw()

CRC<-know2[which(know2$Type=="CRC"),]
HCC<-know2[which(know2$Type=="HCC"),]
NP<-know2[which(know2$Type=="NP"),]
PAAD<-know2[which(know2$Type=="PAAD"),]

CRC$Type<-NULL
HCC$Type<-NULL
NP$Type<-NULL
PAAD$Type<-NULL

densityCRC<-as.data.frame(rowSums(CRC))
densityCRC$tag<-"CRC"
colnames(densityCRC)<-c("Totalcount","Tag")

densityPAAD<-as.data.frame(rowSums(PAAD))
densityPAAD$tag<-"PAAD"
colnames(densityPAAD)<-c("Totalcount","Tag")

densityNP<-as.data.frame(rowSums(NP))
densityNP$tag<-"NP"
colnames(densityNP)<-c("Totalcount","Tag")

densityHCC<-as.data.frame(rowSums(HCC))
densityHCC$tag<-"HCC"
colnames(densityHCC)<-c("Totalcount","Tag")

density<-rbind(densityCRC,densityPAAD,densityNP,densityHCC)

##################################################
###Abundance of circRNA in exosomes
library("ggplot2")
ggplot(data=density,aes(x=Totalcount,fill=factor(Tag),alpha=0.4))+
geom_density(stat="density")+
scale_x_continuous(name="Total count of circRNA in individual samples",expand=c(0,0))+
scale_y_continuous(name="Density",expand=c(0,0))+
theme_classic()

cleanCRC<-as.data.frame(table(colSums(as.matrix(CRC)>0)))
cleanHCC<-as.data.frame(table(colSums(as.matrix(HCC)>0)))
cleanNP<-as.data.frame(table(colSums(as.matrix(NP)>0)))
cleanPAAD<-as.data.frame(table(colSums(as.matrix(PAAD)>0)))

cleanCRC<-cleanCRC[-1,]
cleanHCC<-cleanHCC[-1,]
cleanNP<-cleanNP[-1,]
cleanPAAD<-cleanPAAD[-1,]

cleanCRC$Var1<-as.character(cleanCRC$Var1)
cleanHCC$Var1<-as.character(cleanHCC$Var1)
cleanNP$Var1<-as.character(cleanNP$Var1)
cleanPAAD$Var1<-as.character(cleanPAAD$Var1)
cleanCRC$Var1<-as.numeric(cleanCRC$Var1)
cleanHCC$Var1<-as.numeric(cleanHCC$Var1)
cleanNP$Var1<-as.numeric(cleanNP$Var1)
cleanPAAD$Var1<-as.numeric(cleanPAAD$Var1)

cleanCRC$tag<-"CRC"
cleanHCC$tag<-"HCC"
cleanNP$tag<-"NP"
cleanPAAD$tag<-"PAAD"
cleanCRC$Var1p<-cleanCRC$Var1/12
cleanHCC$Var1p<-cleanHCC$Var1/21
cleanNP$Var1p<-cleanNP$Var1/32
cleanPAAD$Var1p<-cleanPAAD$Var1/14

cleanall<-rbind(cleanCRC,cleanHCC,cleanNP,cleanPAAD)
cleanall$logfreq<-log10(cleanall$Freq)

###Sparsity of circRNA in exosomes
ggplot(data=cleanall, mapping=aes(x=Var1p, y=logfreq,group=tag)) +
  geom_line(aes(color=tag))+
  scale_x_continuous(name="Proportion of samples haboring the circRNA")+
  scale_y_continuous(name="log10(Types of the circRNA)")+
  theme_classic() 

################################################
###expression of pan-cancer or pan-normal tissue-enriched circRNA
allcirclist<-read.csv("allcirclist.csv",header=F)

select<-intersect(allcirclist$V1,colnames(know2))

length(select)

select<-c(select,"Type")

know21<-know2[,select]

my_comparisons<-list(c("NP","PAAD"),c("NP","CRC"),c("NP","HCC"))

for(i in 1:(ncol(know21)-1)){
  df1<-know21[,c(i,ncol(know21))]
  name<-colnames(df1)[1]
  p<-ggboxplot(df1, x="Type", y=name, fill = "Type",order=c("NP","PAAD","CRC","HCC"))+
 stat_compare_means(comparisons = my_comparisons, label = "p.signif")+#label这里表示选择显著性标记（星号） 
    stat_compare_means(label.y = 50)
p
  myfilename <- paste(name,".pdf",sep="")
  ggsave(filename = myfilename,width = 3, height = 5)
}

knowmut<-know2[,colnames(know2)[grep("chr1_224952669",colnames(know2))]]

library("ggpubr")
p<-ggboxplot(know2, x="Type", y="chr1_224952669_224955098_+", fill = "Type",order=c("NP","PAAD","CRC","HCC"))+
 stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
    stat_compare_means(label.y = 50)
p

##################################################
###abundance of cancer-specific, pan-cancer tissue-enriched and pan-normal tissue-enriched circRNA

know21t<-as.data.frame(t(know21))

know21t<-know21t[allcirclist$V1,]

know21t<-cbind(know21t,allcirclist$label)

know21t<-know21t[grep("chr",rownames(know21t)),]

know21cs<-know21t[which(know21t$'allcirclist$label'=="cancer-specific"),]
know21ce<-know21t[which(know21t$'allcirclist$label'=="cancer-enriched"),]
know21ne<-know21t[which(know21t$'allcirclist$label'=="normal-enriched"),]

know21cs$'allcirclist$label'<-NULL
know21ce$'allcirclist$label'<-NULL
know21ne$'allcirclist$label'<-NULL

for(i in 1:79){
know21cs[,i]<-as.numeric(know21cs[,i])
}

for(i in 1:79){
know21ce[,i]<-as.numeric(know21ce[,i])
}

for(i in 1:79){
know21ne[,i]<-as.numeric(know21ne[,i])
}

csrowmean<-as.data.frame(rowMeans(know21cs))
cerowmean<-as.data.frame(rowMeans(know21ce))
nerowmean<-as.data.frame(rowMeans(know21ne))

csrowmean

cscolmean<-as.data.frame(colMeans(know21cs))
colnames(cscolmean)<-"totalcount"
cscolmean$tag<-"cancer-specific"
cecolmean<-as.data.frame(colMeans(know21ce))
colnames(cecolmean)<-"totalcount"
cecolmean$tag<-"cancer-enriched"
necolmean<-as.data.frame(colMeans(know21ne))
colnames(necolmean)<-"totalcount"
necolmean$tag<-"normal-enriched"

cscolmean$sample<-know21$Type
cecolmean$sample<-know21$Type
necolmean$sample<-know21$Type

density<-rbind(cscolmean,cecolmean,necolmean)

densityCRC<-density[which(density$sample=="CRC"),]
densityPAAD<-density[which(density$sample=="PAAD"),]
densityNP<-density[which(density$sample=="NP"),]
densityHCC<-density[which(density$sample=="HCC"),]

densityCRC$totalcount<-log10(densityCRC$totalcount+1)
densityPAAD$totalcount<-log10(densityPAAD$totalcount+1)
densityNP$totalcount<-log10(densityNP$totalcount+1)
densityHCC$totalcount<-log10(densityHCC$totalcount+1)

library(ggplot2)
library(lattice)
library(Rmisc)

library("ggplot2")
p1<-ggplot(data=densityNP,aes(x=totalcount,fill=factor(tag),alpha=0.4))+
geom_density(stat="density")+
scale_x_continuous(name="log10(Total count of circRNA in NP samples+1)",expand=c(0,0))+
scale_y_continuous(name="Density",expand=c(0,0))+
theme_classic()
p2<-ggplot(data=densityPAAD,aes(x=totalcount,fill=factor(tag),alpha=0.4))+
geom_density(stat="density")+
scale_x_continuous(name="log10(Total count of circRNA in PAAD samples+1)",expand=c(0,0))+
scale_y_continuous(name="Density",expand=c(0,0))+
theme_classic()
p3<-ggplot(data=densityCRC,aes(x=totalcount,fill=factor(tag),alpha=0.4))+
geom_density(stat="density")+
scale_x_continuous(name="log10(Total count of circRNA in CRC samples+1)",expand=c(0,0))+
scale_y_continuous(name="Density",expand=c(0,0))+
theme_classic()
p4<-ggplot(data=densityHCC,aes(x=totalcount,fill=factor(tag),alpha=0.4))+
geom_density(stat="density")+
scale_x_continuous(name="log10(Total count of circRNA in HCC samples+1)",expand=c(0,0))+
scale_y_continuous(name="Density",expand=c(0,0))+
theme_classic()

multiplot(p1,p2,p3,p4,cols=2)

