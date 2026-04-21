#Paper UME  Figure3

library("data.table")
library("ggplot2")
library("ggpmisc")
library("patchwork")
library("data.table")
library("ggridges")
library("ggplotify")
library("cowplot")
library("gridBase")





Sys.setenv(LANG = "en")
setwd("L:/MAVE/Experiments/Hydra_coverage/OL")
library("colourvalues")

#file='MATRIX.max.testsubplusGATK.V7.hmp_IBS.dist.txt'

file='MATRIX.max.testsubplusGATK.V7.2.hmp_IBS.dist.txt'

#colrf='Color_key_Phred_cutoff5.txt'
colrf='Color_key_Phred_cutoff6.txt'


cuan=1



distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
#copy  row names and create colnames as it is a matrix
colnames(distancesAll) <- rownames(distancesAll)

#read colors file ##header is: Taxon	MDSGroup	ColByMDSGroup
col <- read.table(colrf,header=T,sep="\t",as.is=T, fill=TRUE)

colrs<-colour_values(1:13, palette="rainbow", alpha=255)



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

outerClass  <- "group1"
innerClass <- "group2"

legend.inner <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),innerClass],levels = unique(col[,innerClass])))
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,"color",sep="_")],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]


# plot solution


  main<-""
  # svg(out)
  i<-1
  j <- i+1
  x <- fit$points[,i]
  y <- fit$points[,j]

  par(mfrow = c(1, 1), mar=c(3,3,0.1,0.7),bg=NA)
  plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=3,bg=paste(color.outer,"FF",sep=""),pch=21,cex=1.3,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*0),cex.axis=1,cex.lab=1, mgp = c(1.5, 0.5, 0))
  legend(0.11,0.01 ,legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=0.75,pt.cex=1,bty="n",y.intersp=2)
  pmds1<-recordPlot()
  #dev.off()








######################### distros

#!/usr/bin/env Rscript
Sys.setenv(LANG = "en")
#create a MDS plot from a matrix file

#setwd("L:/MAVE/phred")
setwd("L:/MAVE/SCRATCH_BKP_2021/V7")

colrf='Color_key_Phred_cutoff2.txt'
cuan=1

vals=0
algs=0
cuts=0


par(mfrow = c(1, 1), mar=c(4,5,2,1))

testsa<-list("A"="HM", "B"=0, "C"=2, "D"=4, "E"=5, "F"=6, "G"=7,"H"=8,"I"=10, "J"=1, "K"=3, "M"=9)
for (cur in names(testsa)) {
  qv<- testsa[[cur]]
  cutoff=qv
  #print(cur)
  #print(cutoff)


  tests<-list("A"="min" ,"B"="max", "C"="median", "D"="mean")
  #tests<-list("A"="max")


  for (cur in names(tests)) {

    alg<- tests[[cur]]
    #print(cur)
    #print(alg)
    file<-paste("MATRIX.",alg,".testsub10HMdB.V7.IBS.txt",sep="")
    #### load distance matrix
    distancesAll <- read.table(file,sep="\t",header=T,row.names=1,as.is=T)
    #copy  row names and create colnames as it is a matrix
    colnames(distancesAll) <- rownames(distancesAll)

    #read colors file ##header is: Taxon	MDSGroup	ColByMDSGroup
    col <- read.table(colrf,header=T,sep="\t",as.is=T,comment.char = "@")

    #col <- subset(col, col$Group2==cutoff | col$Group2=="HM"  )
    col <- subset(col, col$Group2==cutoff  )
    col <- subset(col, col$Group1!="Ky21")

    # Things that can be done to filter the key file and choose some groups

    #### update list of values to plot based on intersection between color list and individuals in file
    distances <- distancesAll[which(rownames(distancesAll)%in%col$Taxon==T),which(rownames(distancesAll)%in%col$Taxon==T)]
    col <- col[which(col$Taxon%in%rownames(distances)),]
    ##### this helps to scale the plot based on the maximum values
    ds <- as.matrix(distances)
    #plot results in MDS space
    fit <- cmdscale(ds,k=2,eig=T) # k is the number of dim
    #fit # view results

    i=1
    j <- i+1
    k <-j+1
    l <-k+1

    z<-round((fit$eig[i]/sum(fit$eig))*100,1)
    w<-round((fit$eig[j]/sum(fit$eig))*100,1)
    r<-round((fit$eig[k]/sum(fit$eig))*100,1)
    m<-round((fit$eig[l]/sum(fit$eig))*100,1)

    val=z+w+r+m

    vals <- c(vals, val)
    algs <- c(algs, alg)
    cuts <- c(cuts, cutoff)


  }
}


df<- data.frame(vals,algs,cuts)

library(ggplot2)
library(hrbrthemes)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 5
colrs = gg_color_hue(n)


#install.packages("hrbrthemes")

dfx=subset(df,!vals==0)


legend_title=" "
p2<-ggplot(dfx, aes(x=cuts, y=vals, color=algs, group=algs)) +
  scale_color_manual(values=colrs)+
  labs(color = " ")+
  geom_line(linetype="aa",aes(color=algs),size=0.6)+
  geom_point(size=1.8,shape=1,stroke = 1)+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))+
  theme(legend.key = element_rect(colour = NA, fill = NA))+
  theme(axis.title=element_text(size=12))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(legend.position = c(0.2, 0.9))+
  theme(legend.background = element_rect(fill='transparent'))+
  labs(y="Variance explained (%)", x = "Normalized Phred cutoff")+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_x_discrete(limits = c("HM","0", "1","2","3","4","5","6","7","8","9","10"))
p2




library("colourvalues")
library("data.table")

setwd("L:/MAVE/SCRATCH_BKP_2021/V7")






par(mfrow = c(1, 1), mar=c(4,5,2,1))

##################Calculations for plot J
dHM <- nrow(na.omit(fread("hmp321_agpv5_chr10.unimputed.vcf.gz.282.qual.lst")))
tests<-list("A"=0 ,"B"=0.1, "C"=0.5, "D"=1, "E"=2,"EA"=3,"F"=4, "G"=5, "H"=6, "I"=7, "J"=8, "JA"=9, "K"=10)
column.names <- c("max", "min", "mean", "median","sum","cutoff")
row.names <- c(0,0.1,0.5,1,2,3,4,5,6,7,8,9,10)
E<-matrix(data=0, nrow=13, ncol=6, dimnames=list(row.names, column.names));

counter=0
for (cur in names(tests)) {
  counter = counter + 1
  qv<- tests[[cur]]
  print(cur)
  print(qv)

  phred=qv

  dmax <- nrow(na.omit(fread(paste("TESTsubset.poli.10.Union.max.db.quals.",phred,".282.lst",sep=""))))
  dmin <- nrow(na.omit(fread(paste("TESTsubset.poli.10.Union.min.db.quals.",phred,".282.lst",sep=""))))
  dmean <- nrow(na.omit(fread(paste("TESTsubset.poli.10.Union.mean.db.quals.",phred,".282.lst",sep=""))))
  dmedian <- nrow(na.omit(fread(paste("TESTsubset.poli.10.Union.median.db.quals.",phred,".282.lst",sep=""))))
  dsum <- nrow(na.omit(fread(paste("TESTsubset.poli.10.Union.sum.db.quals.",phred,".282.lst",sep=""))))

  E[counter,1]=dmax
  E[counter,2]=dmin
  E[counter,3]=dmean
  E[counter,4]=dmedian
  E[counter,5]=dsum
  E[counter,6]=phred

}

#E.df <- as.data.frame(t(E))
E.df <- as.data.frame(E)

###calculate max and min values of the distributions

mxvy<-max(E.df$max,E.df$min,E.df$mean,E.df$median,E.df$sum)
mxvx<-max(E.df$max,E.df$min,E.df$mean,E.df$median,E.df$sum)
mnvy<-min(E.df$max,E.df$min,E.df$mean,E.df$median,E.df$sum)
mnvx<-min(E.df$max,E.df$min,E.df$mean,E.df$median,E.df$sum)
a<-as.numeric(mxvy)
mxvy<-a




#pdf("NSNPS_V7.pdf",width=8,height=8)



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 5
colrs = gg_color_hue(n)

#colrs<-colour_values(1:5, palette="plasma", alpha=255)


par(mfrow = c(1, 1), mar=c(3,3,1,0.5))

plot(y=E.df$max,x=E.df$cutoff, type= "b",col=colrs[1],lwd=1.5,xlab="Normalized Phred cutoff", ylab="Number of variants", main="",cex.axis=1, cex.lab=1,cex.main=1.5,mgp = c(1.5, 0.5, 0))
lines(y=E.df$mean,x=E.df$cutoff, type= "b",col=colrs[2],lwd=1.5)
lines(y=E.df$median,x=E.df$cutoff, type= "b",col=colrs[3],lwd=1.5)
lines(y=E.df$min,x=E.df$cutoff, type= "b",col=colrs[4],lwd=1.5)
legend(7.5,1400000, legend=c("max","mean","median","min","HM3"), col=colrs, pt.cex=1, pch=1,cex=1, bty="n", y.intersp=2)
abline(h=dHM, col=colrs[5],lwd = 2,lty=2)
p3<-recordPlot()

dHM


#### COW PLOTING

setwd("L:/MAVE/Figures")

FIGB<-plot_grid(p2, p3, labels = c('B', 'C'), label_size = 20, hjust = 0, vjust = 1, scale= 0.97)
FIGB

FIGF<-plot_grid(pmds1,FIGB,label_size = 20,labels = c('A',''), hjust = 0, vjust = 1, scale= 1,ncol=1)
FIGF

pdf("Vallebueno2022_Fig_3.pdf",width=6.69,height = 6.5)
FIGF
dev.off()


