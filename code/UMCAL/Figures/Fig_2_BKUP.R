#Paper UME  Figure2

library("data.table")
library("ggplot2")
library("ggpmisc")
library("patchwork")
#install.packages("ggridges")
library("data.table")
library("ggridges")


#!/usr/bin/env Rscript
Sys.setenv(LANG = "en")
#create a MDS plot from a matrix file

setwd("L:/MAVE/Experiments/Hydra_coverage/INT1")


file='MERGE_merges.consensus.1.geno.hmp_IBS.dist.txt'
colrf='Color_key_con1_NC.lst'

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

col<-col[order(col$group1),]

legend.inner <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),innerClass],levels = unique(col[,innerClass])))
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.inner.legend <- c(col[match(legend.inner,as.vector(col[,innerClass])),paste(innerClass,"color",sep="_")],rep("#FFFFFF",length(legend.outer)))
color.outer.legend <- c(rep("#FFFFFF",length(legend.inner)),col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.inner,legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]


##plot

out=gsub(".txt",paste(".mds.",cuan,".svg",sep=""), file)
tti=gsub(".geno.hmp_IBS.dist.txt"," callers", file)
tti=gsub("Merge_Hydratest_","", tti)
tti=gsub(".geno.consensus."," ", tti)
main <-  paste("MDS", tti)
main<-""
# svg(out)
i=1
j <- i+1
x <- fit$points[,i]
y <- fit$points[,j]

library("ggplotify")
library("cowplot")
library("gridBase")

setwd("L:/MAVE/SCRATCH_BKP_2021/V7")
main=""
#par(mfrow = c(1, 1), mar=c(5,5,0.1,0.1))

#svg("MDS_GATK_bias_DP.svg",width=12,height=10)


text=legend
par(mfrow = c(1, 1), mar=c(3,3,0.1,0.1),bg=NA)
plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=3,bg=paste(color.outer,"FF",sep=""),pch=21,cex=1.3,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*0),cex.axis=1,cex.lab=1, mgp = c(1.5, 0.5, 0))
legend(x = "bottomleft",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=0.75,pt.cex=1,bty="n",y.intersp=1.6)
#dev.off()



pmds<-recordPlot()




######################### distros

setwd("L:/MAVE/phred")
dataM <- fread("FIN_SUBSET_1.ALL.Merge.db.phredquant.norm")
dataF <- subset(dataM, dataM$qual<310)
datax <- subset(dataF, dataF$caller!="")
dataF=datax
p2<-ggplot(dataF, aes(x = dataF$qual, y = dataF$caller)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size=12))+
  theme(legend.position="none")+
  #theme(legend.title = element_text(size=30))+
  #theme(legend.text = element_text(size=20))+
  #theme(legend.position = 'bottom')+
  geom_density_ridges(aes(fill =  dataF$caller)) +
  labs(y=" ", x = "Phred",fill = " ")+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#18c4ff", "#a7bc86"))
p3<-ggplot(dataF, aes(x = dataF$norm, y = dataF$caller)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size=12))+
  theme(legend.position="none")+
  #theme(legend.title = element_text(size=30))+
  #theme(legend.text = element_text(size=20))+
  #theme(legend.position = 'bottom')+
  labs(y= " ", x = "Norm.pval",fill=" ")+
  geom_density_ridges(aes(fill =  dataF$caller)) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#18c4ff", "#a7bc86"))












#### COW PLOTING

setwd("L:/MAVE/Figures")

FIGB<-plot_grid(p2, p3, labels = c('B', 'C'), label_size = 20)
FIGB

FIGF<-plot_grid(pmds,FIGB,label_size = 20,labels = c('A',''), hjust = 0, vjust = 1, scale= 1,ncol=1)
FIGF



pdf("Vallebueno2022_Fig_2.pdf",width=6.69,height = 6.5)
FIGF
dev.off()
