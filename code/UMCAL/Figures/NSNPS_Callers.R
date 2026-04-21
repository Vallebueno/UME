setwd("L:/MAVE/SCRATCH_BKP_2021/V7")
library("data.table")
library("colourvalues")
#FILE="RES_DIST2cont.lst"
FILE="RES_DIST2contV7.lst"
#https://www.tenderisthebyte.com/blog/2019/04/25/rotating-axis-labels-in-r/
#https://stackoverflow.com/questions/14604439/plot-multiple-boxplot-in-one-graph



colrs<-colour_values(1:11, palette="matlab_like", alpha=255)


nnm=c("","","","","","","max","","","","","","","","","","","mean","","","","","","","","","","","median","","","","","","","","","","","min","","","","","")

data=fread(FILE)
#dat <- subset(data, data$value<0.1)
dat <- subset(data, data$name!= "sum")

dats <- subset(dat, dat$group == 10)


par(mfrow = c(1, 1), mar=c(4,5,2,1))
out="IBSDist2control.svg"
svg(out,width=15,height=10)

#out="IBSDist2control.pdf"
#pdf(out,width=15,height=10)




boxplot(dat$value ~ dat$group + dat$name, data = dat, col=colrs,boxwex=1.1,xlab="", ylab="IBS distance",cex.axis=1.5,cex.lab=1.5,
        xaxs = FALSE, xaxt = "n"
)

axis(side = 1, at = seq_along(nnm), labels = nnm, tick = FALSE, cex.axis=1.6)
legend("topright", fill = colrs, legend = c(0,1,2,3,4,5,6,7,8,9,10), horiz = F, title="N.Phred cutoff",cex=1.3)
stripchart(dat$value ~ dat$group + dat$name, method="jitter", pch=19, col =1, vertical = TRUE, add=TRUE)


dev.off()