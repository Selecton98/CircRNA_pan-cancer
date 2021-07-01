
xcell<-read.table(file="xCell_mioncocircgeneexp_xCell_0029062021.txt",sep="\t",header=T)

rownames(xcell)<-xcell[,1]

xcell<-xcell[,-1]

allmeta<-read.csv(file="mioncometadata.csv")

allmeta<-allmeta[,-c(1,2)]

rownames(allmeta)<-allmeta$ID

colnames(xcell)<-gsub("\\.","-",colnames(xcell))

colnames(xcell)<-gsub("-csv","\\.csv",colnames(xcell))

see<-intersect(colnames(xcell),allmeta$ID)

allmeta2<-allmeta[colnames(xcell),]

load(file="allmioncocirc.Rda") #too big for github.available on e-mail request.

rownames(all)<-all$circRNA
all$circRNA<-NULL

all<-all[,allmeta2$ID]

allmeta2can<-allmeta2[which(allmeta2$Type=="cancer"),]
allmeta2nor<-allmeta2[which(allmeta2$Type=="normal"),]

anno<-as.data.frame(allmeta2[,2])

rownames(anno)<-allmeta2$ID
colnames(anno)<-"anno"

#circRNA:all
#microenviroment:xcell
#metadata:allmeta2;allmeta2can;allmeta2nor

xcell<-xcell[c(23,37,17,62,64,4,5,10,18,19,20,21,31,35,39,41,42,44,47),]

library(pheatmap)
pheatmap(xcell,annotation_col=anno,cluster_rows = T,cluster_col=T, scale="row",color = colorRampPalette(c("navy","white","purple"))(100),cellwidth=2,cellheight=5,fontsize = 5,border=F)

enriched<-read.csv(file="allcirclist.csv",header=T)

enriched<-enriched[which(enriched$label!="cancer-specific"),]

library(tidyr)
enriched<-separate(enriched,V1,into= c("chr","start","end","strand"),sep= "\\_")

enriched$circRNA<-paste(enriched$chr,enriched$start,enriched$end,sep="_")

enriched$circRNA

all<-all[enriched$circRNA,]

xcellcor<-as.data.frame(t(xcell))
circcor<-as.data.frame(t(all))
circcor<-circcor[rownames(xcellcor),]

xcellcan<-xcell[,allmeta2can$ID]
allcan<-all[,allmeta2can$ID]

xcellcancor<-as.data.frame(t(xcellcan))
circcancor<-as.data.frame(t(allcan))
circcancor<-circcor[rownames(xcellcancor),]

xcellnor<-xcell[,allmeta2nor$ID]
allnor<-all[,allmeta2nor$ID]

xcellnorcor<-as.data.frame(t(xcellnor))
circnorcor<-as.data.frame(t(allnor))
circnorcor<-circcor[rownames(xcellnorcor),]

circcor<-circcor[,which(colSums(circcor)>0),]
circcancor<-circcancor[,which(colSums(circcancor)>0),]
circnorcor<-circnorcor[,which(colSums(circnorcor)>0),]

anno<-enriched[,c(6,7,8)]
rownames(anno)<-anno$circRNA
anno$circRNA<-NULL

anno$V2<-NULL

anno$label<-gsub("-",".",anno$label)

ann_colors=list(label=c(cancer.enriched="tomato",normal.enriched="turquoise"))

table(anno$label)


xcellcancor2<-xcellcancor[,c(23,37,17,62,64,4,5,10,18,19,20,21,31,35,39,41,42,44,47)]

library(psych)
res<-corr.test(xcellcancor, circcancor, use = "pairwise",method="spearman",adjust="bonferroni", alpha=.05)
pheatmap(res$r,cluster_rows = T,annotation_col=anno,annotation_colors=ann_colors,cluster_col=T,display_numbers = matrix(ifelse(res$p <= 0.01, "**", ifelse(res$p <= 0.05 ,"*"," ")), nrow(res$p)), color = colorRampPalette(c("steelblue","white","firebrick"))(100),cellwidth=10,cellheight=10,fontsize =10,border=F)

##transcription-driven circRNA
load(file="mioncocircgeneexp.Rda")

geneexp2<-as.data.frame(t(geneexp))
geneexp3<-geneexp2[rownames(circcancor),]
geneexp4<-geneexp2[rownames(circnorcor),]

bindcan<-cbind(circcancor,geneexp3)
bindnor<-cbind(circcancor,geneexp4)

library(ggplot2)
library(ggpubr)
library(ggExtra)

p1 <- ggplot(bindcan, aes(DNAH14,chr1_224952669_224968874)) + 
  geom_point(color="#6baed6") +
  geom_smooth(method=lm , color = "#756bb1", fill = "#cbc9e2" ,se=FALSE) +
  stat_cor(data=bindcan, method = "spearman")+
  geom_hline(yintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  theme_classic()
             
p1<-ggMarginal(p1,type="histogram",fill="#756bb1")

p2 <- ggplot(bindnor, aes(DNAH14,chr1_224952669_224968874)) + 
  geom_point(color="#6baed6") +
  geom_smooth(method=lm , color = "#756bb1", fill = "#cbc9e2" ,se=FALSE) +
  stat_cor(data=bindnor, method = "spearman")+
  geom_hline(yintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  theme_classic()
             
p2<-ggMarginal(p2,type="histogram",fill="turquoise")


p3 <- ggplot(bindcan, aes(DNAH14,chr1_224952669_224974153)) + 
  geom_point(color="#6baed6") +
  geom_smooth(method=lm , color = "#756bb1", fill = "#cbc9e2" ,se=FALSE) +
  stat_cor(data=bindcan, method = "spearman")+
  geom_hline(yintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  theme_classic()
             
p3<-ggMarginal(p3,type="histogram",fill="#756bb1")

p4 <- ggplot(bindnor, aes(DNAH14,chr1_224952669_224974153)) + 
  geom_point(color="#6baed6") +
  geom_smooth(method=lm , color = "#756bb1", fill = "#cbc9e2" ,se=FALSE) +
  stat_cor(data=bindnor, method = "spearman")+
  geom_hline(yintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  theme_classic()
             
p4<-ggMarginal(p4,type="histogram",fill="turquoise")
