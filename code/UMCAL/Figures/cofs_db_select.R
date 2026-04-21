#!/usr/bin/env Rscript

library("colourvalues")
Sys.setenv(LANG = "en")

library("data.table")


colrs<-colour_values(1:7, palette="matlab_like", alpha=255)



#create a MDS plot from a matrix file

#setwd("L:/MAVE/phred")
setwd("L:/MAVE/Experiments/DB_COF_alldata")


#colrf='Colors_ND_capvsWGS_6_PTC.txt'
colrf='Colors_ALL_August.txt'

ordls<-'ord_leg.lst'


df<-read.table(ordls, header=F, sep= "\t")
ordl<-df$V1
ordl



tests<-list("A"="cof5" ,"B"="cof6", "C"="cof8", "D"="cof10")

par(mfrow = c(2, 2),mar=c(4,5,2,1))
for (cur in names(tests)) {

  alg<- tests[[cur]]
  print(cur)
  print(alg)
  file<-paste(alg,"_Chr1_tlone_Mergedb.vcf_IBS.dist.txt",sep="")


  #file='cof5_Chr1_tlone_Mergedb.vcf_IBS.dist.txt'

  distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
  colnames(distancesAll) <- rownames(distancesAll)
  col <- read.table(colrf,header=T,sep="\t",as.is=T)

  distances <- distancesAll[which(rownames(distancesAll)%in%col$Taxon==T),which(rownames(distancesAll)%in%col$Taxon==T)]
  col <- col[which(col$Taxon%in%rownames(distances)),]
  dsx <- as.matrix(distances)
  ds<-dsx[complete.cases(dsx), complete.cases(dsx)]


  fit <- cmdscale(ds,k=3,eig=T) # k is the number of dim

  outerClass  <- "Site"
  innerClass <- "Type"


  legend.inner <- unique(col[match(rownames(fit$points[]),col$Taxon,),innerClass])
  legend.inner <- ordl[which(ordl%in%legend.inner==T)]
  legend.outer <- unique(col[match(rownames(fit$points[]),col$Taxon,),outerClass])
  legend.outer <- ordl[which(ordl%in%legend.outer==T)]

  color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste("color",sep="_",innerClass)],rep("#FFFFFF",length(legend.outer)))
  color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste("color",sep="_",outerClass)])
  legend <- c(legend.inner,legend.outer)
  color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste("color",sep="_",innerClass)]
  color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste("color",sep="_",outerClass)]


  for (i in 1:1) {
    j <- i+1
    x <- fit$points[,i]
    y <- fit$points[,j]
    plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=alg,col=color.inner,lwd=4,bg=paste(color.outer,"",sep=""),pch=21,cex=2.7,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1.6, cex.lab=1.7,cex.main=2)
    legend(x = "topright",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=1,pt.cex=2)

  }
}
#dev.off()











#####################################################


file='cof5_Chr1_tlone_Mergedb.vcf_IBS.dist.txt'

distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
colnames(distancesAll) <- rownames(distancesAll)
col <- read.table(colrf,header=T,sep="\t",as.is=T)

distances <- distancesAll[which(rownames(distancesAll)%in%col$Taxon==T),which(rownames(distancesAll)%in%col$Taxon==T)]
col <- col[which(col$Taxon%in%rownames(distances)),]
dsx <- as.matrix(distances)
ds<-dsx[complete.cases(dsx), complete.cases(dsx)]


fit <- cmdscale(ds,k=2,eig=T) # k is the number of dim

outerClass  <- "Site"
innerClass <- "Type"

legend.inner <- unique(col[match(rownames(fit$points[]),col$Taxon,),innerClass])
legend.inner <- ordl[which(ordl%in%legend.inner==T)]
legend.outer <- unique(col[match(rownames(fit$points[]),col$Taxon,),outerClass])
legend.outer <- ordl[which(ordl%in%legend.outer==T)]

color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste("color",sep="_",innerClass)],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste("color",sep="_",outerClass)])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste("color",sep="_",innerClass)]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste("color",sep="_",outerClass)]


for (i in 1:1) {
  j <- i+1
  x <- fit$points[,i]
  y <- fit$points[,j]
  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.inner,lwd=4,bg=paste(color.outer,"",sep=""),pch=21,cex=2.7,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1.6, cex.lab=1.7,cex.main=2)
  legend(x = "topright",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=1,pt.cex=2)
}
