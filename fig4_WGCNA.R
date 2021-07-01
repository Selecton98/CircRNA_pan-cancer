load("IDCSCmatrix.Rda")
load("IDCSCmeta.Rda")
rownames(matrix3)<-matrix3$circRNA
matrix3$circRNA<-NULL
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.6)
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.7)
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.8)
quantile(rowSums(as.matrix(matrix3[1:584])> 0),0.9)
matrix3q<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=1&rowSums(as.matrix(matrix3)> 0)<2),]
matrix3q6<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=2&rowSums(as.matrix(matrix3)> 0)<4),]
matrix3q7<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=4&rowSums(as.matrix(matrix3)> 0)<7),]
matrix3q8<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=7&rowSums(as.matrix(matrix3)> 0)<20),]
matrix3q9<-matrix3[which(rowSums(as.matrix(matrix3)> 0)>=20),]
dat<-matrix3q9 #according to which quantile to be investigated
############################################################
###WGCNA co-expression analysis
library(WGCNA)
library(reshape2)
library(stringr)
dataExpr<-as.data.frame(t(dat))
csv<-colnames(dataExpr)
type = "unsigned"
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average") 
pdf("hcluster.pdf")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
abline(h = 5000, col = "red")
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,networkType="unsigned",verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
pdf("power.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()
power = sft$powerEstimate
power
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type = "unsigned",9, 18),
               ifelse(nSamples<30, ifelse (type = "unsigned",8, 16),
                        ifelse(nSamples<40, ifelse(type = "unsigned",7, 14),
                               ifelse(type = "unsigned",6, 12))       
                 )
  )
}
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = "unsigned", minModuleSize = 5, 
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = "pearson", 
                       loadTOMs=TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
save(net,file="net.Rda")
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
pdf("dendro.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
dim(TOM)
save(probes,file="probes.Rda")
save(moduleColors,file="moduleColors.Rda")
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf("eigengene.pdf")
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
rownames<-rownames(dataExpr)
rownames<-as.data.frame(rownames)
colnames(rownames)<-"Samples"
color1<-merge(rownames,meta4,all=F)
cancer<-gsub("cancer","1",color1$Type)  
cancer<-gsub("normal","0",cancer)  
normal<-gsub("cancer","0",color1$Type)  
normal<-gsub("normal","1",normal)  
datTraits<-cbind(rownames,cancer,normal)
colnames(datTraits)<-c("sample","cancer","normal")
rownames(datTraits)<-datTraits[,1]
datTraits<-datTraits[,-1]
modTraitCor = cor(MEs_col, datTraits, use = "p")
modTraitP = corPvalueStudent(modTraitCor,nSamples)
modTraitCor = cor(MEs_col, datTraits, use = "p")
modTraitP = corPvalueStudent(modTraitCor,nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf("moduletrait.pdf")
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(datTraits), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.3, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
############################################################
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
