library("colourvalues")
library("data.table")

setwd("L:/MAVE/SCRATCH_BKP_2021/V7")


out="Phred_Dist_V7.svg"
svg(out,width=15,height=15)

#out="Phred_Dist_V7.pdf"
#pdf(out,width=10,height=10)



par(mfrow = c(1, 1), mar=c(4,5,2,1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 5
colrs = gg_color_hue(n)
phred=5
dmax <- na.omit(fread(paste("TESTsubset.poli.10.Union.max.db.quals.",phred,".282.lst",sep="")))
dmin <- na.omit(fread(paste("TESTsubset.poli.10.Union.min.db.quals.",phred,".282.lst",sep="")))
dmean <- na.omit(fread(paste("TESTsubset.poli.10.Union.mean.db.quals.",phred,".282.lst",sep="")))
dmedian <- na.omit(fread(paste("TESTsubset.poli.10.Union.median.db.quals.",phred,".282.lst",sep="")))
#dsum <- na.omit(fread(paste("TESTsubset.poli.10.Union.sum.db.quals.",phred,".282.lst",sep="")))

###calculate max and min values of the distributions
mxvy<-max(density(dmax$V1)$y,density(dmin$V1)$y,density(dmean$V1)$y,density(dmedian$V1)$y)#,density(dsum$V1)$y)
mxvx<-max(density(dmax$V1)$x,density(dmin$V1)$x,density(dmean$V1)$x,density(dmedian$V1)$x)#,density(dsum$V1)$x)
#mnvy<-min(density(dmax$V1)$y,density(dmin$V1)$y,density(dmean$V1)$y,density(dmedian$V1)$y)#,density(dsum$V1)$y)
#mnvx<-min(density(dmax$V1)$x,density(dmin$V1)$x,density(dmean$V1)$x,density(dmedian$V1)$x)#,density(dsum$V1)$x)


plot(density(dmax$V1),col=colrs[1],lwd=2,xlab="Phred", ylab="Density", main="Distribution of phred by algoritm",cex.axis=1.2, cex.lab=1.4,cex.main=1.5,ylim=c(0,mxvy))
lines(density(dmin$V1),col=colrs[4],lwd=2)
lines(density(dmean$V1),col=colrs[2],lwd=2)
lines(density(dmedian$V1),col=colrs[3],lwd=2)
legend("topleft", legend=c("max","mean","median", "min"), col=colrs[1:5], pt.cex=2.5, pch=15,cex=1.5)
#mtext(text="I", xpd=NA, side=2, adj=1, font=2,cex=2.2,line=2,las=1,padj=0,at=mxvy*1.2)






dev.off()
