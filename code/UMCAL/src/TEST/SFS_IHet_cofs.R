
library("colourvalues")
library("data.table")

setwd("L:/MAVE/TEST")

library("colourvalues")
colrs<-colour_values(1:11, palette="viridis", alpha=255)


par(mfrow = c(1, 1), mar=c(4,5,2,1))



#tests<-list("A"="min" ,"B"="max", "C"="median", "D"="mean")


#for (cur in names(tests)) {
#  alg<- tests[[cur]]
#  print(cur)
#  print(alg)

  dHM <- na.omit(fread("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.1.mins.cuantex"))

  plot(density(dHM$V1),col=colrs[1],lwd=2,xlab="minor AF", ylab="Density", main="", cex.axis=1.2, cex.lab=1.4,ylim=c(0,500),xlim=c(0,0.5),cex.main=1.5,,cex.axis=1.6, cex.lab=1.7,cex.main=2)
  #title(alg, line = 0.6, cex.main = 2)


  testsq<-list("2"=2,"3"=4 ,"4"=8 ,"5"=16, "6"=32, "7"=64, "8"=96, "9"=112, "10"=120, "11"=124)

  for (curq in names(testsq)) {
    cof<- testsq[[curq]]
    print(curq)
    print(cof)
    n=as.numeric(curq)

   # df <- na.omit(fread(paste("TESTsubset.poli.10.Union.",alg,".db.",phred,".282.SFS",sep="")))
    df <- na.omit(fread(paste("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.",cof,".mins.cuantex",sep="")))

    lines(density(df$V1),col=colrs[n],lwd=2)
  }
  legend("topright", legend=c("1","2","4","8","16","32","64", "96", "112", "120", "124"), col=colrs, pt.cex=3, pch=15,cex=1.7)



#}

dev.off()



##########SFS per grp

colrs<-colour_values(1:7, palette="viridis", alpha=255)
par(mfrow = c(1, 1), mar=c(4,5,2,1))


dHM <- na.omit(fread("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.I1.grp.cuantex"))

plot(density(dHM$V1),col=colrs[1],lwd=2,xlab="minor AF", ylab="Density", main="", cex.axis=1.2, cex.lab=1.4,ylim=c(0,500),xlim=c(0,0.5),cex.main=1.5,,cex.axis=1.6, cex.lab=1.7,cex.main=2)
#title(alg, line = 0.6, cex.main = 2)


testsq<-list("2"="I2","3"="C" ,"4"="A1" ,"5"="A2", "6"="T", "7"="P")

for (curq in names(testsq)) {
  cof<- testsq[[curq]]
  print(curq)
  print(cof)
  n=as.numeric(curq)

  # df <- na.omit(fread(paste("TESTsubset.poli.10.Union.",alg,".db.",phred,".282.SFS",sep="")))
  df <- na.omit(fread(paste("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.",cof,".grp.cuantex",sep="")))

  lines(density(df$V1),col=colrs[n],lwd=2)
}
legend("topright", legend=c("I1","I2","C","A1","A2","T","P"), col=colrs, pt.cex=3, pch=15,cex=1.7)







#### SFS of variant comming from grp

colrs<-colour_values(1:7, palette="viridis", alpha=255)
par(mfrow = c(1, 1), mar=c(4,5,2,1))


dHM <- na.omit(fread("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.I1.cgrp.cuantey"))

plot(density(dHM$V1),col=colrs[1],lwd=2,xlab="minor AF", ylab="Density", main="", cex.axis=1.2, cex.lab=1.4,ylim=c(0,150),xlim=c(0,0.5),cex.main=1.5,,cex.axis=1.6, cex.lab=1.7,cex.main=2)
#title(alg, line = 0.6, cex.main = 2)


testsq<-list("2"="I2","3"="C" ,"4"="A1" ,"5"="A2", "6"="T", "7"="P")

for (curq in names(testsq)) {
  cof<- testsq[[curq]]
  print(curq)
  print(cof)
  n=as.numeric(curq)

  # df <- na.omit(fread(paste("TESTsubset.poli.10.Union.",alg,".db.",phred,".282.SFS",sep="")))
  df <- na.omit(fread(paste("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.",cof,".cgrp.cuantey",sep="")))

  lines(density(df$V1),col=colrs[n],lwd=2)
}
legend("topright", legend=c("I1","I2","C","A1","A2","T","P"), col=colrs, pt.cex=3, pch=15,cex=1.7)


####### OBS-EXP

par(mfrow = c(1, 1), mar=c(4,5,2,1))

nms=matrix(rep(0,300),nrow=300)

colrs<-colour_values(1:10, palette="viridis", alpha=255)


cof=2
df <- na.omit(fread(paste("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.",cof,".IHET.cuantey.eo",sep="")))
plot(NULL, xlim=c(0,0.055), ylim=c(0,0.04), xlab="Expected",ylab="Observed" )
#lines(lowess(df$V3,df$V4), col=colrs[1], lwd=4.0)

mod1=lm(df$V3~df$V4)
modsum = summary(mod1)
abline(mod1, col=colrs[1], lwd=4.0)
r2 = modsum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
nms[2,1] = round(r2, digits=3)

testsq<-list("2"=4 ,"3"=8 ,"4"=16, "5"=32, "6"=64, "7"=96, "8"=112, "9"=120, "10"=124)

for (curq in names(testsq)) {
  cof<- testsq[[curq]]
  print(curq)
  print(cof)
  n=as.numeric(curq)

  # df <- na.omit(fread(paste("TESTsubset.poli.10.Union.",alg,".db.",phred,".282.SFS",sep="")))
  df <- na.omit(fread(paste("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.",cof,".IHET.cuantey.eo",sep="")))

  mod1=lm(df$V3~df$V4)
  modsum = summary(mod1)
  abline(mod1, col=colrs[n], lwd=4.0)
  r2 = modsum$adj.r.squared
  nms[cof,1]= round(r2, digits=3)

}


legend("topleft", legend=c(paste("2",nms[2,1],sep=" r2 = "),paste("4",nms[4,1],sep=" r2 = "),paste("8",nms[8,1],sep=" r2 = "),paste("16",nms[16,1],sep=" r2 = "),paste("32",nms[32,1],sep=" r2 = "),paste("64",nms[64,1],sep=" r2 = "),paste("96",nms[96,1],sep=" r2 = "),paste("112",nms[112,1],sep=" r2 = "),paste("120",nms[120,1],sep=" r2 = "),paste("124",nms[124,1],sep=" r2 = ")), col=colrs, pt.cex=3, pch=15,cex=1.2)


par(mfrow = c(2, 5), mar=c(4,5,2,1))


testsq<-list("1"=2,"2"=4 ,"3"=8 ,"4"=16, "5"=32, "6"=64, "7"=96, "8"=112, "9"=120, "10"=124)

for (curq in names(testsq)) {
  cof<- testsq[[curq]]
  print(curq)
  print(cof)
  n=as.numeric(curq)



df <- na.omit(fread(paste("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.",cof,".IHET.cuantey.eo",sep="")))
plot(x=df$V4,y=df$V3, xlim=c(0,0.055), ylim=c(0,0.04), xlab="Expected",ylab="Observed", main=paste("HET Cutoff",cof) )
#lines(lowess(df$V3,df$V4), col=colrs[1], lwd=4.0)
abline(lm(df$V3~df$V4), col="red", lwd=4.0)
}


#####################singletons from each group


colrs<-colour_values(1:7, palette="viridis", alpha=255)
par(mfrow = c(1, 1), mar=c(4,5,2,1))


dHM <- na.omit(fread("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.I1.sgrp.cuantex"))

plot(density(dHM$V1),col=colrs[1],lwd=2,xlab="minor AF", ylab="Density", main="", cex.axis=1.2, cex.lab=1.4,ylim=c(0,150),xlim=c(0,0.5),cex.main=1.5,,cex.axis=1.6, cex.lab=1.7,cex.main=2)
#title(alg, line = 0.6, cex.main = 2)


testsq<-list("2"="I2","3"="C" ,"4"="A1" ,"5"="A2", "6"="T", "7"="P")

for (curq in names(testsq)) {
  cof<- testsq[[curq]]
  print(curq)
  print(cof)
  n=as.numeric(curq)

  # df <- na.omit(fread(paste("TESTsubset.poli.10.Union.",alg,".db.",phred,".282.SFS",sep="")))
  df <- na.omit(fread(paste("1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.",cof,".cgrp.cuantey",sep="")))

  lines(density(df$V1),col=colrs[n],lwd=2)
}
legend("topright", legend=c("I1","I2","C","A1","A2","T","P"), col=colrs, pt.cex=3, pch=15,cex=1.7)


##################### Dist sim het vs general dist


dHM <- na.omit(fread("1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf..DHET.cuantex"))
dHM1 <- na.omit(fread("1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf_IBS.dist.txt.DHET.cuantex"))

testsq<-list("1"=2,"2"=4 ,"3"=8 ,"4"=16, "5"=32, "6"=64, "7"=96, "8"=112)

par(mfrow = c(4, 2), mar=c(4,5,2,1))
for (curq in names(testsq)) {
  cof<- testsq[[curq]]
  print(curq)
  print(cof)
  n=as.numeric(curq)


  xdHM<-dHM[which(dHM$V1%in%(cof))]
  xdHM1<-dHM1[which(dHM1$V1%in%(cof))]

 # xdHM<-dHM[which(dHM$V1%in%(100))]

d<-density(xdHM$V2)
d1<-density(xdHM1$V2)
  cof<-96
  cof<-112
plot(d,ylim=c(0,150),xlim=c(0,0.6),col="red",main=paste(cof,"inds with Hets"),yaxt="n",xaxt="n",ylab="",xlab="")
  axis(1,cex.axis=1.5)
  axis(2,cex.axis=1.5)
  mtext("Genetic distance", side=1, line=2.2, cex=1)
  mtext("Density", side=2, line=2.2, cex=1)
ls
lines(d1,col="blue")
legend("topright", legend=c("Dist_Hets","Dist_Rand"), col=c("red","blue"), pt.cex=2, pch=15,cex=1.7)

}

#hist(dHM$V1, breaks = 100, col=rgb(1,0,0,0.5),freq = FALSE , xlab="Distance" , ylab="Frequency" , main="" )
#hist(dHM1$V1, col=rgb(0,0,1,0.5) ,add=T,freq = FALSE)




####### MDS

library("colourvalues")

colrs<-colour_values(1:14, palette="rainbow", alpha=255)
colrs<-colour_values(1:12, palette="viridis", alpha=255)
colrs<-colour_values(1:5, palette="greys", alpha=255)

colrs<-colour_values(1:8, palette="greys", alpha=255)

colrs<-colour_values(1:32, palette="rainbow", alpha=255)


#/groups/swarts/lab/MAVE/MAIZEDB/ID.key

setwd("L:/MAVE/TEST")

file='1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf_IBS.dist.txt'
#colrf='Color_key_NAM.txt'
#colrf='Color_key_CAPL.txt'
#colrf='Color_key_grpsnoT.txt'
#colrf='Color_key_grpsnoT.txt'
#colrf='Color_key_grpANC.txt'
colrf='Color_key_countr_modern.txt'
#colrf='Color_key_CL_L.txt'
#colrf='Color_key_grpsNoplusancient1.txt'

par(mfrow = c(1, 1), mar=c(4,5,2,1))
cuan=1



distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)




#copy  row names and create colnames as it is a matrix
colnames(distancesAll) <- rownames(distancesAll)

#read colors file ##header is: Taxon	MDSGroup	ColByMDSGroup
col <- read.table(colrf,header=T,sep="\t",as.is=T, fill=TRUE)
colrs<-colour_values(1:8, palette="viridis", alpha=255)



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
fit <- cmdscale(ds,k=2,eig=T) # k is the number of dim
#fit # view results


#MDS Group for color

outerClass  <- "group2"
innerClass <- "group1"

legend.inner <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),innerClass],levels = unique(col[,innerClass])))
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,"color",sep="_")],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]


# plot solution

for (i in 1:1) {
  # svg(out)
  j <- i+1
  x <- fit$points[,i]
  y <- fit$points[,j]

#  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=4,bg=paste(color.outer,"FF",sep=""),pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.outer,lwd=6,bg=color.inner,pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

  legend(x = "topright",legend = legend,col=color.outer.legend,pt.bg=color.inner.legend,pt.lwd = 3,pch=21,cex=1.7,pt.cex=4)
#  legend(x = "topright",legend = legend,col=color.outer.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=2,pt.cex=3)
  #dev.off()
}







###############################  MDS ancient +modern grays



####### MDS

library("colourvalues")

colrs<-colour_values(1:12, palette="rainbow", alpha=255)
colrs<-colour_values(1:12, palette="viridis", alpha=255)
colrs<-colour_values(1:5, palette="greys", alpha=255)

colrs<-colour_values(1:14, palette="greys", alpha=255)

colrs<-colour_values(1:32, palette="rainbow", alpha=255)


#/groups/swarts/lab/MAVE/MAIZEDB/ID.key

setwd("L:/MAVE/TEST")

file='1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf_IBS.dist.txt'
colrf='Color_key_NAM.txt'
#colrf='Color_key_CAPL.txt'
#colrf='Color_key_grpsnoT.txt'
#colrf='Color_key_grpsnoT.txt'
colrf='Color_key_grpANC.txt'
#colrf='Color_key_countr_modern.txt'
#colrf='Color_key_CL_L.txt'
#colrf='Color_key_grpsNoplusancient1.txt'

par(mfrow = c(1, 1), mar=c(4,5,2,1))
cuan=1



distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)




#copy  row names and create colnames as it is a matrix
colnames(distancesAll) <- rownames(distancesAll)

#read colors file ##header is: Taxon	MDSGroup	ColByMDSGroup
col <- read.table(colrf,header=T,sep="\t",as.is=T, fill=TRUE)
colrs<-colour_values(1:8, palette="viridis", alpha=255)



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
innerClass <- "group1"


legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))



color.outer.legend <- c(col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]


# plot solution

#for (i in 1:1) {
  # svg(out)
i <-1
  j <- i+1
  x <- fit$points[,i]
  y <- fit$points[,j]

  #  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=4,bg=paste(color.outer,"FF",sep=""),pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.outer,lwd=6,bg=color.inner,pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

  legend(x = "topright",legend = legend,col=color.outer.legend,pt.lwd = 3,pch=16,cex=1.7,pt.cex=4)
  #  legend(x = "topright",legend = legend,col=color.outer.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=2,pt.cex=3)
  #dev.off()
#}

fit$points[,1]

Xa<-as.data.frame(fit$points)

fit$points


#############suplementary


setwd("L:/MAVE/TEST")

file='1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf_IBS.dist.txt'
colrf='Color_key_NAM.txt'

par(mfrow = c(1, 1), mar=c(4,5,2,1))
cuan=1



distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)




#copy  row names and create colnames as it is a matrix
colnames(distancesAll) <- rownames(distancesAll)

#read colors file ##header is: Taxon	MDSGroup	ColByMDSGroup
col <- read.table(colrf,header=T,sep="\t",as.is=T, fill=TRUE)
colrs<-colour_values(1:8, palette="viridis", alpha=255)



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
fit <- cmdscale(ds,k=2,eig=T) # k is the number of dim
#fit # view results


#MDS Group for color

outerClass  <- "group2"
innerClass <- "group1"

legend.inner <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),innerClass],levels = unique(col[,innerClass])))
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,"color",sep="_")],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]


# plot solution

for (i in 1:1) {
  # svg(out)
  j <- i+1
  x <- fit$points[,i]
  y <- fit$points[,j]

  #  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=4,bg=paste(color.outer,"FF",sep=""),pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.outer,lwd=6,bg=color.inner,pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=2.5, cex.lab=2.5,cex.main=2)

  legend(x = "topright",legend = legend,col=color.outer.legend,pt.bg=color.inner.legend,pt.lwd = 3,pch=21,cex=1.7,pt.cex=4)
  #  legend(x = "topright",legend = legend,col=color.outer.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=2,pt.cex=3)
  #dev.off()
}




###########################sample correlation age ,, coverage etc


library("colourvalues")
library("data.table")

setwd("L:/MAVE/TEST")


colrs<-colour_values(1:11, palette="viridis", alpha=255)

par(mfrow = c(1, 1), mar=c(4,5,2,1))


data <- na.omit(fread("Ancient_samples_info.txt"))


par(mfrow = c(3, 4), mar=c(4,5,2,1))

#### PC1 vs mean coverage

mma<-max(data$PC1)
mmb<-min(data$PC1)
nma<-max(data$mean_genome_coverage)
nmb<-min(data$mean_genome_coverage)

plot(x=data$mean_genome_coverage,y=data$PC1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="mean genome coverage",ylab="PC1", main="PC1 vs mean_genome_coverage")
abline(lm(data$PC1~as.factor(data$Dataset)+as.factor(data$Geography)+data$Age+data$mean_genome_coverage), col="red", lwd=4.0)
xx<-lm(data$PC1~as.factor(data$Dataset)+as.factor(data$Geography)+data$Age+data$mean_genome_coverage)
summary(xx)

mma<-max(data$mean_genome_coverage)
mmb<-min(data$mean_genome_coverage)
nma<-max(data$PC1)
nmb<-min(data$PC1)

plot(x=data$PC1,y=data$mean_genome_coverage, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean genome coverage", main="PC1 vs mean_genome_coverage")
abline(lm(data$PC1~as.factor(data$Dataset)+as.factor(data$Geography)+data$Age+data$mean_genome_coverage), col="red", lwd=4.0)

xx<-lm(data$PC1~as.factor(data$Dataset)+as.factor(data$Geography)+data$Age+data$mean_genome_coverage)

summary(xx)






#### PC1 vs lenght coverage

mma<-max(data$bases_mincov1)
mmb<-min(data$bases_mincov1)
nma<-max(data$PC1)
nmb<-min(data$PC1)

plot(x=data$PC1,y=data$bases_mincov1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean lenght coverage", main="PC1 vs mean_length_coverage")
abline(lm(data$bases_mincov1~data$PC1), col="red", lwd=4.0)

#### PC1 vs mean_genome_depth

mma<-max(data$mean_genome_depth)
mmb<-min(data$mean_genome_depth)
nma<-max(data$PC1)
nmb<-min(data$PC1)

plot(x=data$PC1,y=data$mean_genome_depth, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean lenght coverage", main="PC1 vs mean_genome_depth")
abline(lm(data$mean_genome_depth~data$PC1), col="red", lwd=4.0)


#### PC1 vs age

mma<-max(data$Age)
mmb<-min(data$Age)
nma<-max(data$PC1)
nmb<-min(data$PC1)

plot(x=data$PC1,y=data$Age, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean lenght coverage", main="PC1 vs age")
abline(lm(data$Age~data$PC1), col="red", lwd=4.0)


##############PC2

#### PC2 vs mean coverage

mma<-max(data$mean_genome_coverage)
mmb<-min(data$mean_genome_coverage)
nma<-max(data$PC2)
nmb<-min(data$PC2)

plot(x=data$PC2,y=data$mean_genome_coverage, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean genome coverage", main="PC2 vs mean_genome_coverage")
abline(lm(data$mean_genome_coverage~data$PC2), col="red", lwd=4.0)

#### PC2 vs lenght coverage

mma<-max(data$bases_mincov1)
mmb<-min(data$bases_mincov1)
nma<-max(data$PC2)
nmb<-min(data$PC2)

plot(x=data$PC2,y=data$bases_mincov1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC2",ylab="mean lenght coverage", main="PC2 vs mean_length_coverage")
abline(lm(data$bases_mincov1~data$PC2), col="red", lwd=4.0)

#### PC2 vs mean_genome_depth

mma<-max(data$mean_genome_depth)
mmb<-min(data$mean_genome_depth)
nma<-max(data$PC2)
nmb<-min(data$PC2)

plot(x=data$PC2,y=data$mean_genome_depth, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC2",ylab="mean lenght coverage", main="PC2 vs mean_genome_depth")
abline(lm(data$mean_genome_depth~data$PC2), col="red", lwd=4.0)


#### PC2 vs age

mma<-max(data$Age)
mmb<-min(data$Age)
nma<-max(data$PC2)
nmb<-min(data$PC2)

plot(x=data$PC2,y=data$Age, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC2",ylab="mean lenght coverage", main="PC2 vs age")
abline(lm(data$Age~data$PC2), col="red", lwd=4.0)




##############PC3

#### PC3 vs mean coverage

mma<-max(data$mean_genome_coverage)
mmb<-min(data$mean_genome_coverage)
nma<-max(data$PC3)
nmb<-min(data$PC3)

plot(x=data$PC3,y=data$mean_genome_coverage, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean genome coverage", main="PC3 vs mean_genome_coverage")
abline(lm(data$mean_genome_coverage~data$PC3), col="red", lwd=4.0)

#### PC3 vs lenght coverage

mma<-max(data$bases_mincov1)
mmb<-min(data$bases_mincov1)
nma<-max(data$PC3)
nmb<-min(data$PC3)

plot(x=data$PC3,y=data$bases_mincov1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC3",ylab="mean lenght coverage", main="PC3 vs mean_length_coverage")
abline(lm(data$bases_mincov1~data$PC3), col="red", lwd=4.0)

#### PC3 vs mean_genome_depth

mma<-max(data$mean_genome_depth)
mmb<-min(data$mean_genome_depth)
nma<-max(data$PC3)
nmb<-min(data$PC3)

plot(x=data$PC3,y=data$mean_genome_depth, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC3",ylab="mean lenght coverage", main="PC3 vs mean_genome_depth")
abline(lm(data$mean_genome_depth~data$PC3), col="red", lwd=4.0)


#### PC3 vs age

mma<-max(data$Age)
mmb<-min(data$Age)
nma<-max(data$PC3)
nmb<-min(data$PC3)

plot(x=data$PC3,y=data$Age, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC3",ylab="mean lenght coverage", main="PC3 vs age")
abline(lm(data$Age~data$PC3), col="red", lwd=4.0)


#############Modern


###########################sample correlation age ,, coverage etc


library("colourvalues")
library("data.table")

setwd("L:/MAVE/TEST")


colrs<-colour_values(1:11, palette="viridis", alpha=255)

par(mfrow = c(1, 1), mar=c(4,5,2,1))


data <- na.omit(fread("Modern_samples_info.txt"))


par(mfrow = c(3, 3), mar=c(4,5,2,1))



#### PC1 vs mean coverage

mma<-max(data$mean_genome_coverage)
mmb<-min(data$mean_genome_coverage)
nma<-max(data$PC1)
nmb<-min(data$PC1)

plot(x=data$PC1,y=data$mean_genome_coverage, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean genome coverage", main="PC1 vs mean_genome_coverage")
abline(lm(data$mean_genome_coverage~data$PC1), col="red", lwd=4.0)

#### PC1 vs lenght coverage

mma<-max(data$bases_mincov1)
mmb<-min(data$bases_mincov1)
nma<-max(data$PC1)
nmb<-min(data$PC1)

plot(x=data$PC1,y=data$bases_mincov1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean lenght coverage", main="PC1 vs mean_length_coverage")
abline(lm(data$bases_mincov1~data$PC1), col="red", lwd=4.0)

#### PC1 vs mean_genome_depth

mma<-max(data$mean_genome_depth)
mmb<-min(data$mean_genome_depth)
nma<-max(data$PC1)
nmb<-min(data$PC1)

plot(x=data$PC1,y=data$mean_genome_depth, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean lenght coverage", main="PC1 vs mean_genome_depth")
abline(lm(data$mean_genome_depth~data$PC1), col="red", lwd=4.0)



##############PC2

#### PC2 vs mean coverage

mma<-max(data$mean_genome_coverage)
mmb<-min(data$mean_genome_coverage)
nma<-max(data$PC2)
nmb<-min(data$PC2)

plot(x=data$PC2,y=data$mean_genome_coverage, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean genome coverage", main="PC2 vs mean_genome_coverage")
abline(lm(data$mean_genome_coverage~data$PC2), col="red", lwd=4.0)

#### PC2 vs lenght coverage

mma<-max(data$bases_mincov1)
mmb<-min(data$bases_mincov1)
nma<-max(data$PC2)
nmb<-min(data$PC2)

plot(x=data$PC2,y=data$bases_mincov1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC2",ylab="mean lenght coverage", main="PC2 vs mean_length_coverage")
abline(lm(data$bases_mincov1~data$PC2), col="red", lwd=4.0)

#### PC2 vs mean_genome_depth

mma<-max(data$mean_genome_depth)
mmb<-min(data$mean_genome_depth)
nma<-max(data$PC2)
nmb<-min(data$PC2)

plot(x=data$PC2,y=data$mean_genome_depth, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC2",ylab="mean lenght coverage", main="PC2 vs mean_genome_depth")
abline(lm(data$mean_genome_depth~data$PC2), col="red", lwd=4.0)





##############PC3

#### PC3 vs mean coverage

mma<-max(data$mean_genome_coverage)
mmb<-min(data$mean_genome_coverage)
nma<-max(data$PC3)
nmb<-min(data$PC3)

plot(x=data$PC3,y=data$mean_genome_coverage, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC1",ylab="mean genome coverage", main="PC3 vs mean_genome_coverage")
abline(lm(data$mean_genome_coverage~data$PC3), col="red", lwd=4.0)

#### PC3 vs lenght coverage

mma<-max(data$bases_mincov1)
mmb<-min(data$bases_mincov1)
nma<-max(data$PC3)
nmb<-min(data$PC3)

plot(x=data$PC3,y=data$bases_mincov1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC3",ylab="mean lenght coverage", main="PC3 vs mean_length_coverage")
abline(lm(data$bases_mincov1~data$PC3), col="red", lwd=4.0)

#### PC3 vs mean_genome_depth

mma<-max(data$mean_genome_depth)
mmb<-min(data$mean_genome_depth)
nma<-max(data$PC3)
nmb<-min(data$PC3)

plot(x=data$PC3,y=data$mean_genome_depth, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="PC3",ylab="mean lenght coverage", main="PC3 vs mean_genome_depth")
abline(lm(data$mean_genome_depth~data$PC3), col="red", lwd=4.0)

############################################################# corrected Linear regression


library("colourvalues")
library("data.table")

setwd("L:/MAVE/TEST")


colrs<-colour_values(1:11, palette="viridis", alpha=255)

par(mfrow = c(1, 1), mar=c(4,5,2,1))


data <- na.omit(fread("Ancient_Mod_samples_info.txt"))


par(mfrow = c(3, 4), mar=c(4,5,2,1))

#### PC1 vs mean coverage


mma<-max(data$PC1)
mmb<-min(data$PC1)
nma<-max(data$mean_genome_coverage)
nmb<-min(data$mean_genome_coverage)



xlmo<-lm(data$PC1~as.factor(data$Dataset)+as.factor(data$Geography)+as.factor(data$AM)+as.factor(data$info)+data$Age+data$mean_genome_coverage)
summary(xlmo)


plot(x=data$mean_genome_coverage,y=data$PC1, xlim=c(nmb,nma), ylim=c(mmb,mma), xlab="mean genome coverage",ylab="PC1", main="PC1 vs mean_genome_coverage")
abline(fitted(xlmo), col="red", lwd=4.0)

fitted(xlmo)
