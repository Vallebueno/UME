#Paper UME  Figure5

library("data.table")
library("ggplot2")
library("ggpmisc")
library("patchwork")
#install.packages("ggridges")
library("data.table")
library("ggridges")

library("colourvalues")

setwd("L:/MAVE/TEST")

file='1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf_IBS.dist.txt'
#file='GBSModern_ancient_nocap_filt2.vcf.gz_IBS.dist.txt'
#colrf='Color_key_NAM.txt'
#colrf='Color_key_CAPL.txt'
#colrf='Color_key_grpsnoT.txt'
#colrf='Color_key_grpsnoT.txt'
#colrf='Color_key_grpANC.txt'
colrf='Color_ANcient_VIP.txt'
#colrf='Color_key_countr_modern.txt'
#colrf='Color_key_CL_L.txt'
#colrf='Color_key_grpsNoplusancient1.txt'
#colrf='head_nocap2.txt'
par(mfrow = c(1, 1), mar=c(4,5,2,1))
cuan=1



distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)


#copy  row names and create colnames as it is a matrix
colnames(distancesAll) <- rownames(distancesAll)

#read colors file ##header is: Taxon	MDSGroup	ColByMDSGroup
col <- read.table(colrf,header=T,sep="\t",as.is=T, fill=TRUE)


#colrs<-colour_values(1:8, palette="rainbow", alpha=255)
#colrs<-colour_values(1:60, palette="greys", alpha=255)
colrs<-colour_values(1:13, palette="blue2red", alpha=255)


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



names(which(colSums(is.na(ds))>0))
names(which(rowSums(is.na(ds))>0))



#plot results in MDS space
fit <- cmdscale(ds,k=4,eig=T) # k is the number of dim
#fit # view results


#MDS Group for color

outerClass  <- "group2"
innerClass <- "group2"


legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.outer.legend <- c(col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]




#legend.inner <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),innerClass],levels = unique(col[,innerClass])))
#legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
#color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,"color",sep="_")],rep("#FFFFFF",length(legend.outer)))
#color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
#legend <- c(legend.inner,legend.outer)
#color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
#color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]








# plot solution

#for (i in 1:1) {
# svg(out)
i <-1
j <- i+1
x <- fit$points[,i]
y <- fit$points[,j]

#  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=4,bg=paste(color.outer,"FF",sep=""),pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.outer,lwd=6,bg=color.inner,pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

legend(x = "topright",legend = legend,col=color.outer.legend,pt.lwd = 3,pch=16,cex=0.8,pt.cex=2)



#legend(x = "topleft",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=0.75,pt.cex=1,bty="n",y.intersp=2)






#  legend(x = "topright",legend = legend,col=color.outer.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=2,pt.cex=3)
#dev.off()
#}





fit$points[,1]





############################

install.packages("data.table")
install.packages("ggplot2")
library(data.table)
setwd("L:/MAVE/TEST")

dataM <- fread("coverage_inGBS.txt")

dfa <- subset(dataM, dataM$Group == "Ancient")
dfm <- subset(dataM, dataM$Group == "Modern")


hist(dfa$`Proportion Missing`, xlim=c(0,1), ylim=c(0,130),col=rgb(0,0,1,0.5),xlab="Proportion missing", ylab="Frequency",main="",cex.axis=1.2,cex.lab=1.4)
hist(dfm$`Proportion Missing`, breaks=30, xlim=c(0,300), col=rgb(1,0,0,0.5), add=T)


legend("topright", legend=c("Ancient","Modern"), col=c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), pt.cex=4, pch=15,cex=2)


