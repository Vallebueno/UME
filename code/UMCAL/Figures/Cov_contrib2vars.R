setwd("L:/MAVE/Experiments/Hydra_coverage/OL")
library("data.table")
FILE="files.mergein3.lst.tlone.Merge.db.7.poli.qual0.Union.V7-1.max.db.qualst.txt"

data=fread(FILE, header=T)



hist(data$`30`,breaks=100, main="", xlab = 'Phred', col=rgb(1,0,0,0.5) )
hist(data$`5`, breaks=100,add=T , col=rgb(0,0,1,0.5) )
hist(data$`3`, breaks=100,add=T , col=rgb(0,1,1,0.5) )
hist(data$`0.5`, breaks=100, add=T , col=rgb(0,1,0,0.5) )

legend("topright", legend=c("30X","5X", "3X", "0.5X" ), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(0,1,1,0.5), rgb(0,1,0,0.5)), pt.cex=2, pch=15 )

abline(v=5,col="red", lwd=3, lty=2)


