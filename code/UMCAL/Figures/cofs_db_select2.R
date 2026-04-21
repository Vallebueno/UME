#!/usr/bin/env Rscript

library("colourvalues")
Sys.setenv(LANG = "en")

library("data.table")


colrs<-colour_values(1:7, palette="matlab_like", alpha=255)



#create a MDS plot from a matrix file

setwd("L:/MAVE/Experiments/DB_COF_alldata")


colrf='Colors_Nteo_August.txt'


plot.new()
tests<-list("1"="cof1" ,"2"="cof2", "3"="cof3", "4"="cof4", "5"="cof5", "6"="cof6", "7"="cof7", "8"="cof8", "9"="cof9", "10"="cof10", "20"="cof20")

par(mfrow = c(4, 3),mar=c(4,5,2,1))
for (cur in names(tests)) {

  alg<- tests[[cur]]
  print(cur)
  print(alg)
 # file<-paste(alg,"_Chr1_tlone_Mergedb.vcf_IBS.dist.txt",sep="")
  file<-paste("files_10.",cur,"_r2_out_IBS.dist.txt",sep="")
  
  #file='cof5_Chr1_tlone_Mergedb.vcf_IBS.dist.txt'

  distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
  colnames(distancesAll) <- rownames(distancesAll)
  col <- read.table(colrf,header=T,sep="\t",as.is=T)

  distances <- distancesAll[which(rownames(distancesAll)%in%col$Taxon==T),which(rownames(distancesAll)%in%col$Taxon==T)]
  col <- col[which(col$Taxon%in%rownames(distances)),]
  ds <- as.matrix(distances)



  fit <- cmdscale(ds,k=2,eig=T) # k is the number of dim

  
  outerClass  <- "group2"
  innerClass <- "group2"


  legend.inner <- unique(col[match(rownames(fit$points[]),col$Taxon,),innerClass])

  legend.outer <- unique(col[match(rownames(fit$points[]),col$Taxon,),outerClass])

  color.outer.legend <- c(col[match(legend,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
  
  

  for (i in 1:1) {
    j <- i+1
    x <- fit$points[,i]
    y <- fit$points[,j]
   # plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=alg,col=color.inner,lwd=4c,pch=21,cex=2.7,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1.6, cex.lab=1.7,cex.main=1)
    #legend(x = "topright",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=1,pt.cex=2)

   
   plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.outer,lwd=2,bg=color.inner,pch=21,cex=0.7,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*0.02),cex.axis=1.2, cex.lab=1.2,cex.main=0.1)
   
    
    
  }
}
#dev.off()
plot.new()
legend("left",legend = legend,col=color.outer.legend, ncol=2,pt.lwd = 0,pch=16,cex=1,pt.cex=1.2,bty="n",y.intersp=1,x.intersp=1)









#####################################################

cur<-1

setwd("L:/MAVE/Experiments/DB_COF_alldata")


colrf='Colors_ALL_August.txt'
file<-paste("files_10.",cur,"_r2_out_IBS.dist.txt",sep="")

#file='cof5_Chr1_tlone_Mergedb.vcf_IBS.dist.txt'

distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
colnames(distancesAll) <- rownames(distancesAll)
col <- read.table(colrf,header=T,sep="\t",as.is=T)

distances <- distancesAll[which(rownames(distancesAll)%in%col$Taxon==T),which(rownames(distancesAll)%in%col$Taxon==T)]
col <- col[which(col$Taxon%in%rownames(distances)),]
#dsx <- as.matrix(distances)
#ds<-dsx[complete.cases(dsx), complete.cases(dsx)]
ds <- as.matrix(distances)

fit <- cmdscale(ds,k=2,eig=T) # k is the number of dim

outerClass  <- "group1"
innerClass <- "group2"


legend.inner <- unique(col[match(rownames(fit$points[]),col$Taxon,),innerClass])

legend.outer <- unique(col[match(rownames(fit$points[]),col$Taxon,),outerClass])


color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,sep="_","color")],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,sep="_","color")])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,sep="_","color")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,sep="_","color")]


for (i in 1:1) {
  j <- i+1
  x <- fit$points[,i]
  y <- fit$points[,j]
  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=alg,col=color.inner,lwd=4,bg=paste(color.outer,"",sep=""),pch=21,cex=2.7,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1.6, cex.lab=1.7,cex.main=2)
  legend(x = "topright",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=1,pt.cex=2)
} 