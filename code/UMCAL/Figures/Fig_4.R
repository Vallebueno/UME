#Paper UME  Figure4

library("data.table")
library("ggplot2")
library("ggpmisc")
library("patchwork")
#install.packages("ggridges")
library("data.table")
library("ggridges")


setwd("L:/MAVE/Experiments/Production")

colrf='heads.key'
file='mpileup.lst.production.db.vcf.hmp_IBS.dist.txt'
cuan=1
distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
#copy  row names and create colnames as it is a matrix
colnames(distancesAll) <- rownames(distancesAll)

#read colors file ##header is: Taxon	MDSGroup	ColByMDSGroup
col <- read.table(colrf,header=T,sep="\t",as.is=T)

# Things that can be done to filter the key file and choose some groups
#col <- col[c(which(col$MDSGroup%in%c("NAM"))),] #"Parviglumis","Mexicana", Remove teosinte,"mexicana","parviglumis" ##Remove groups if desired from col file

#print data about the inputs into out
print("File used for MDS is:")
print(file)
print("Unique values in key file:")
unique(col$MDSGroup)
print("Key colors metrics:")
str(col)

#### update list of values to plot based on intersection between color list and individuals in file
distances <- distancesAll[which(rownames(distancesAll)%in%col$Taxon==T),which(rownames(distancesAll)%in%col$Taxon==T)]
col <- col[which(col$Taxon%in%rownames(distances)),]

##### this helps to scale the plot based on the maximum values
ds <- as.matrix(distances)
#plot results in MDS space
fit <- cmdscale(ds,k=2,eig=T) # k is the number of dim
#fit # view results
#MDS Group for color
#innerClass <- "coverage"
#outerClass <- "group1"
outerClass  <- "group1"
innerClass <- "group2"
# innerClass <- "genotype"
col<-col[order(col$group1),]
legend.inner <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),innerClass],levels = unique(col[,innerClass])))
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,"color",sep="_")],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]

##plot
out=gsub(".txt",paste(".mds.",cuan,".svg",sep=""), file)
tti=gsub(".geno.hmp_IBS.dist.txt"," callers", file)
tti=gsub("Merge_Hydratest_","", tti)
tti=gsub(".geno.consensus."," ", tti)
main <-  paste("MDS", tti)
main<-""
# svg(out)
i=1
j <- i+1
x <- fit$points[,i]
y <- fit$points[,j]
library("ggplotify")
library("cowplot")
library("gridBase")
setwd("L:/MAVE/SCRATCH_BKP_2021/V7")
main=""
text=legend
par(mfrow = c(1, 1), mar=c(3,3,0.1,0.1),bg=NA)
plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=3,bg=paste(color.outer,"FF",sep=""),pch=21,cex=1.3,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*0),cex.axis=1,cex.lab=1, mgp = c(1.5, 0.5, 0))
legend(x = "topleft",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=0.75,pt.cex=1,bty="n",y.intersp=2)
#dev.off()
pmds<-recordPlot()

############ HIST


library("data.table")
library("colourvalues")
setwd("L:/MAVE/Experiments/Hydra_coverage/OL")
dataM <- as.data.frame(fread("MATRIX.max.testsubplusGATK.V7.hmp_SiteSummary.txt"))
dmax<-dataM$`Proportion Missing`
setwd("L:/MAVE/Experiments/Production")
dataP <- as.data.frame(fread("mpileup.lst.production.db.vcf.hmp_SiteSummary.txt"))
dmax<-dataM$`Proportion Missing`
dpro<-dataP$`Proportion Missing`
par(mfrow = c(1, 1), mar=c(3,3,1,0.5))
hist(dpro, breaks = 30, main= "", xlab = 'site proportion missing', col=rgb(1,0,0,0.5),mgp = c(1.5, 0.5, 0),cex.axis=1, cex.lab=1 )
hist(dmax, breaks = 30, add=T , col=rgb(0,0,1,0.5) )
legend(0.65,440000, legend=c("Production","Discovery"), col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)), pt.cex=1, pch=15,cex=1 ,bty="n", y.intersp=2)
p2<-recordPlot()




###########Boxplot

setwd("L:/MAVE/Experiments/Production")
dataP <- as.data.frame(fread("Comp_plot.lst"))
dataP$cov
colrs<-c(rgb(1,0,0,0.5),rgb(0,0,1,0.5))
mis<-dataP$`Proportion Missing`
tax<-dataP$`Taxa Name`
Grop<-dataP$Group
dataP$cov <- factor(dataP$cov, levels=c("0.5x","3x","5x","30x"))

#boxplot(dataP$`Proportion Missing` ~ dataP$Group + dataP$cov +dataP$names , data = dataP, col=colrs,boxwex=1.1,xlab="",ylab="Proportion missing",cex.axis=1.5,cex.lab=1.5,xaxs = FALSE, xaxt = "n")
#boxplot(dataP$`Proportion Missing` ~ dataP$Group + dataP$cov , data = dataP, col=colrs,boxwex=1.1,xlab="",ylab="Proportion missing",cex.axis=1.5,cex.lab=1.5)
par(mfrow = c(1, 1), mar=c(3,3,1,0.5))
boxplot(dataP$`Proportion Missing` ~ dataP$Group + dataP$cov , data = dataP, col=colrs,boxwex=1,xlab="",ylab="",cex.axis=1,cex.lab=1, xaxs = FALSE, xaxt = "n",axes=F)
#boxplot(dataP$`Proportion Missing` ~ dataP$Group + dataP$cov , data = dataP, col=colrs,boxwex=1,xlab="",ylab="",cex.axis=1,cex.lab=1)
#nnm=c("0.5x","","30x","","3x","","5x","")
#nnm=c("0.5x","","","30x","","3x","","5x","")
nnm=c("0.5x","","3x","","5x","","30x","")
axis(side = 1, at = seq_along(nnm), labels = nnm, tick = FALSE, cex.axis=1, line=-1)
axis(side=2,line=0,padj=1, hadj=0.5)
title(xlab = "Taxon coverage", line = 0)            # Add x-axis text
title(ylab = "Proportion missing", line = 0)
legend(5,1, fill = colrs, legend = c("Production","Discovery"), horiz = F, title="",cex=1,pt.cex=1,  y.intersp=2,bty="n")
p3<-recordPlot()



#############hists


setwd("L:/MAVE/Figures")

FIGB<-plot_grid(p2, p3, labels = c('B', 'C'), label_size = 20)
FIGB

FIGF<-plot_grid(pmds,FIGB,label_size = 20,labels = c('A',''), hjust = 0, vjust = 1, scale= 1,ncol=1)
FIGF



pdf("Vallebueno2022_Fig_4.pdf",width=6.69,height = 6.5)
FIGF
dev.off()
