
setwd("L:/MAVE/")
dataP <- as.data.frame(fread("cuant_hmp_int_un.lst"))
dataP$cov
colrs<-c(rgb(1,0,0,0.5),rgb(0,0,1,0.5))

mis<-dataP$num
tax<-dataP$ind
Grop<-dataP$alg
dataP$cov <- factor(dataP$cov, levels=c("0.5x","3x","30x"))



par(mfrow = c(1, 1), mar=c(3,3,1,0.5))
boxplot(dataP$num ~ dataP$alg + dataP$cov , data = dataP, col=colrs,boxwex=1,xlab="",ylab="",cex.axis=1,cex.lab=1, xaxs = FALSE, xaxt = "n",axes=F)
boxplot(dataP$num ~ dataP$alg + dataP$cov , data = dataP, col=colrs,boxwex=1,xlab="",ylab="",cex.axis=1,cex.lab=1 )
title(xlab = "Taxon coverage", line = 2)            # Add x-axis text
title(ylab = "Number of sites", line = 2)
legend(x="bottomright", fill = colrs, legend = c("Intersection","Union"), horiz = F, title="",cex=1,pt.cex=1,  y.intersp=2,bty="n")

aa<-subset(dataP, alg=='Int')

bb<-subset( aa, cov == '0.5x' )

median(bb$num)
