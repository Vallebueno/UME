#!/usr/bin/env Rscript
Sys.setenv(LANG = "en")
setwd("L:/MAVE/Experiments/Hydra_coverage/INT1")

file='MERGE_merges.consensus.1.geno.hmp_IBS.dist.txt'
colrf='Color_key_con1.lst'

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
outerClass  <- "coverage"
innerClass <- "group1"
# innerClass <- "genotype"

legend.inner <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),innerClass],levels = unique(col[,innerClass])))
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,"color",sep="_")],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]


# plot solution

for (i in 1:1) {

  out=gsub(".txt",paste(".mds.",cuan,".svg",sep=""), file)
  tti=gsub(".geno.hmp_IBS.dist.txt"," callers", file)
  tti=gsub("Merge_Hydratest_","", tti)
  tti=gsub(".geno.consensus."," ", tti)
  main <-  paste("MDS", tti)
  main<-""

  #svg("MDS_UNION.svg",width=12,height=10)
 # svg("MDS_INT.svg",width=12,height=10)
  j <- i+1
  x <- fit$points[,i]
  y <- fit$points[,j]

  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=4,bg=paste(color.outer,"FF",sep=""),pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1.6,cex.lab=1.6)
  #text(x, y, labels=rownames(fit$points), xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.acc,pch=19,cex=0.9,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2))
  #text(x, y, labels=rownames(fit$points),main=main,col=color.acc,pch=19,cex=.6,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2))
  legend(x = "bottomright",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=1.25,pt.cex=3)

 # dev.off()
}