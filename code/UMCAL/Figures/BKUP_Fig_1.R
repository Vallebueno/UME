

library("data.table")
library("ggplot2")
library("ggpmisc")
library("patchwork")

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}




setwd("L:/MAVE/SCRATCH_BKP_2021/V7")

#dataU <- fread("TESTsubset.poli.10.Union.max.db.0.phred2dp.lst")

dataU <- na.omit(fread("1.ALL.Merge.10.poli.qual0.Union.V7-1.max.db.0.phred2dp.lst"))
dataU<-subset(dataU, dataU$V2<100)
summary(dataU)
#hist(dataU$V2)
namalg="UMCAL"



setwd("G:/My Drive/Manuscripts/Maize_SNP_DB/Figures/Dist_cov_allsamp")
dataX <- fread("dist_cov_all_samples.txt")

px<-ggplot(dataX, aes(x=dataX$mean_genome_depth)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 35))+
  theme(axis.text.y = element_text(size = 35))+
  theme(axis.title.x = element_text(size=35))+
  theme(axis.title.y = element_text(size=35))+
  theme(legend.title = element_text(size=0))+
  theme(legend.text = element_text(size=35))+
  labs(y="Count", x = "median sample depth",fill = " ")+
  geom_histogram(color="black", fill="blue")

px

setwd("L:/MAVE/phred")

dataB <- fread("freebayes.lst2")
dataG <- fread("gatk.lst2")
dataS <- fread("samtools.lst2")
dataD <- fread("deepvar.lst2")

trsp=0.009 #transparency of dots

ordl=c(namalg,"Deepvar","Samtools","GATK","Freebayes")
################################## p1
samp=10000
phred=20000000000000
mindp=0 ##DP
maxdp=10  ##DP

dataBa <- subset(dataB, dataB$V2<phred) #phred
dataBa <- subset(dataBa, dataBa$V3<maxdp) #dp
dataBa <- subset(dataBa, dataBa$V3>mindp) #dp

dataGa <- subset(dataG, dataG$V2<phred) #phred
dataGa <- subset(dataGa, dataGa$V3<maxdp) #dp
dataGa <- subset(dataGa, dataGa$V3>mindp) #dp

dataSa <- subset(dataS, dataS$V2<phred) #phred
dataSa <- subset(dataSa, dataSa$V3<maxdp) #dp
dataSa <- subset(dataSa, dataSa$V3>mindp) #dp

dataDa <- subset(dataD, dataD$V2<phred) #phred
dataDa <- subset(dataDa, dataDa$V3<maxdp) #dp
dataDa <- subset(dataDa, dataDa$V3>mindp) #dp

dataUa <- subset(dataU, dataU$V1<phred) #phred
dataUa <- subset(dataUa, dataUa$V2<maxdp) #dp
dataUa <- subset(dataUa, dataUa$V2>mindp) #dp


dataBa <- head(dataBa,samp)
dataGa <- head(dataGa,samp)
dataSa <- head(dataSa,samp)
dataDa <- head(dataDa,samp)
dataUa <-head(dataUa,samp)

varu<-rep(namalg,each=samp)
dataUa<-cbind(varu,dataUa)
colnames(dataUa) <- c("V1","V2","V3")

#model=lm(dataDa$V3~dataDa$V2)
#pval <- round(lmp(model),digits=6)
#pval


df_list=list(dataUa,dataDa,dataGa,dataSa,dataBa)
dataF=Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

dataF$V1=factor(dataF$V1,levels=ordl)
dataF=dataF[order(dataF$V1), ]



my.formula <- y ~ x
###lm(phred~dp) is phred predicted by depth





p1<-  ggplot(dataF, aes(x=dataF$V3, y=dataF$V2, color=dataF$V1)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 35))+
  theme(axis.text.y = element_text(size = 35))+
  theme(axis.title.x = element_text(size=35))+
  theme(axis.title.y = element_text(size=35))+
  theme(legend.title = element_text(size=0))+
  theme(legend.text = element_text(size=35))+
  theme(legend.position = c(0.85, 0.85))+
  labs(y="Phred", x = "DP",fill = " ")+
geom_point(aes(color = dataF$V1), alpha=trsp)+
  geom_smooth(method = "lm", fill = NA, se=FALSE,formula = my.formula,size=5)+
    stat_poly_eq(formula = my.formula,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE, size=12)
p1

################################## p2
samp=10000
phred=20000000000000
mindp=30 ##DP
maxdp=100  ##DP

dataBa <- subset(dataB, dataB$V2<phred) #phred
dataBa <- subset(dataBa, dataBa$V3<maxdp) #dp
dataBa <- subset(dataBa, dataBa$V3>mindp) #dp

dataGa <- subset(dataG, dataG$V2<phred) #phred
dataGa <- subset(dataGa, dataGa$V3<maxdp) #dp
dataGa <- subset(dataGa, dataGa$V3>mindp) #dp

dataSa <- subset(dataS, dataS$V2<phred) #phred
dataSa <- subset(dataSa, dataSa$V3<maxdp) #dp
dataSa <- subset(dataSa, dataSa$V3>mindp) #dp

dataDa <- subset(dataD, dataD$V2<phred) #phred
dataDa <- subset(dataDa, dataDa$V3<maxdp) #dp
dataDa <- subset(dataDa, dataDa$V3>mindp) #dp

dataUa <- subset(dataU, dataU$V1<phred) #phred
dataUa <- subset(dataUa, dataUa$V2<maxdp) #dp
dataUa <- subset(dataUa, dataUa$V2>mindp) #dp


dataBa <- head(dataBa,samp)
dataGa <- head(dataGa,samp)
dataSa <- head(dataSa,samp)
dataDa <- head(dataDa,samp)
dataUa <-head(dataUa,samp)

varu<-rep(namalg,each=samp)
dataUa<-cbind(varu,dataUa)


#model=lm(dataDa$V3~dataDa$V2)
#pval <- round(lmp(model),digits=6)
#pval

colnames(dataUa) <- c("V1","V2","V3")
df_list=list(dataUa,dataDa,dataGa,dataSa,dataBa)
dataF2=Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
dataF2$V1=factor(dataF2$V1,levels=ordl)
dataF2=dataF2[order(dataF2$V1), ]
my.formula <- y ~ x
###lm(phred~dp) is phred predicted by depth





p2<-  ggplot(dataF2, aes(x=dataF2$V3, y=dataF2$V2, color=dataF2$V1)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 35))+
  theme(axis.text.y = element_text(size = 35))+
  theme(axis.title.x = element_text(size=35))+
  theme(axis.title.y = element_text(size=35))+
  theme(legend.title = element_text(size=0))+
  theme(legend.text = element_text(size=35))+
  theme(legend.position = c(0.85, 0.85))+
  labs(y="Phred", x = "DP",fill = " ")+
  geom_point(aes(color = dataF2$V1), alpha=trsp)+
  geom_smooth(method = "lm", fill = NA, se=FALSE,formula = my.formula,size=5)+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE, size=12)


p2

################################## p3
samp=10000
phred=20000000000000
mindp=10 ##DP
maxdp=30  ##DP

dataBa <- subset(dataB, dataB$V2<phred) #phred
dataBa <- subset(dataBa, dataBa$V3<maxdp) #dp
dataBa <- subset(dataBa, dataBa$V3>mindp) #dp

dataGa <- subset(dataG, dataG$V2<phred) #phred
dataGa <- subset(dataGa, dataGa$V3<maxdp) #dp
dataGa <- subset(dataGa, dataGa$V3>mindp) #dp

dataSa <- subset(dataS, dataS$V2<phred) #phred
dataSa <- subset(dataSa, dataSa$V3<maxdp) #dp
dataSa <- subset(dataSa, dataSa$V3>mindp) #dp

dataDa <- subset(dataD, dataD$V2<phred) #phred
dataDa <- subset(dataDa, dataDa$V3<maxdp) #dp
dataDa <- subset(dataDa, dataDa$V3>mindp) #dp

dataUa <- subset(dataU, dataU$V1<phred) #phred
dataUa <- subset(dataUa, dataUa$V2<maxdp) #dp
dataUa <- subset(dataUa, dataUa$V2>mindp) #dp


dataBa <- head(dataBa,samp)
dataGa <- head(dataGa,samp)
dataSa <- head(dataSa,samp)
dataDa <- head(dataDa,samp)
dataUa <-head(dataUa,samp)

varu<-rep(namalg,each=samp)
dataUa<-cbind(varu,dataUa)


#model=lm(dataDa$V3~dataDa$V2)
#pval <- round(lmp(model),digits=6)
#pval
colnames(dataUa) <- c("V1","V2","V3")

df_list=list(dataUa,dataDa,dataGa,dataSa,dataBa)
dataF3=Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
dataF3$V1=factor(dataF3$V1,levels=ordl)
dataF3=dataF3[order(dataF3$V1), ]

my.formula <- y ~ x
###lm(phred~dp) is phred predicted by depth



p3<-  ggplot(dataF3, aes(x=dataF3$V3, y=dataF3$V2, color=dataF3$V1)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(legend.title = element_text(size=0))+
  theme(legend.text = element_text(size=20))+
  labs(y="Phred", x = "DP",fill = " ")+
  geom_point(aes(color = dataF3$V1), alpha=trsp)+
  geom_smooth(method = "lm", fill = NA, se=FALSE,formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE)



################################## p4
samp=10000
phred=20000000000000
mindp=1 ##DP
maxdp=70 ##DP

dataBa <- subset(dataB, dataB$V2<phred) #phred
dataBa <- subset(dataBa, dataBa$V3<maxdp) #dp
dataBa <- subset(dataBa, dataBa$V3>mindp) #dp

dataGa <- subset(dataG, dataG$V2<phred) #phred
dataGa <- subset(dataGa, dataGa$V3<maxdp) #dp
dataGa <- subset(dataGa, dataGa$V3>mindp) #dp

dataSa <- subset(dataS, dataS$V2<phred) #phred
dataSa <- subset(dataSa, dataSa$V3<maxdp) #dp
dataSa <- subset(dataSa, dataSa$V3>mindp) #dp

dataDa <- subset(dataD, dataD$V2<phred) #phred
dataDa <- subset(dataDa, dataDa$V3<maxdp) #dp
dataDa <- subset(dataDa, dataDa$V3>mindp) #dp

dataUa <- subset(dataU, dataU$V1<phred) #phred
dataUa <- subset(dataUa, dataUa$V2<maxdp) #dp
dataUa <- subset(dataUa, dataUa$V2>mindp) #dp


dataBa <- head(dataBa,samp)
dataGa <- head(dataGa,samp)
dataSa <- head(dataSa,samp)
dataDa <- head(dataDa,samp)
dataUa <-head(dataUa,samp)

varu<-rep(namalg,each=samp)
dataUa<-cbind(varu,dataUa)
colnames(dataUa) <- c("V1","V2","V3")

#model=lm(dataDa$V3~dataDa$V2)
#pval <- round(lmp(model),digits=6)
#pval


df_list=list(dataUa,dataDa,dataGa,dataSa,dataBa)
dataF4=Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
dataF4$V1=factor(dataF4$V1,levels=ordl)
dataF4=dataF4[order(dataF4$V1), ]
my.formula <- y ~ x
###lm(phred~dp) is phred predicted by depth





p4<-  ggplot(dataF4, aes(x=dataF4$V3, y=dataF4$V2, color=dataF4$V1)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(legend.title = element_text(size=0))+
  theme(legend.text = element_text(size=20))+
  labs(y="Phred", x = "DP",fill = " ")+
  geom_point(aes(color = dataF4$V1), alpha=trsp)+
  geom_smooth(method = "lm", fill = NA, se=FALSE,formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE)




############################ end of plots phred by dist


#!/usr/bin/env Rscript
Sys.setenv(LANG = "en")
#create a MDS plot from a matrix file

setwd("L:/MAVE/Experiments/Hydra_coverage/GATK")


file='Merge_MODERN_GATK_ALLCOV.geno.consensus.1.geno.hmp_IBS.dist.txt'
colrf='Colors_CW4.lst'

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
par(mfrow = c(1, 1), mar=c(5,5,0.1,0.1))

svg("MDS_GATK_bias_DP.svg",width=12,height=10)
plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=4,bg=paste(color.outer,"FF",sep=""),pch=21,cex=3.4,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1.6,cex.lab=1.6)
legend(x = "bottomleft",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=1.55,pt.cex=3)
dev.off()




pmds<-recordPlot()

pdf("MDS_GATK_bias_DP.pdf",width=5,height = 4)
plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main=main,col=color.inner,lwd=3,bg=paste(color.outer,"FF",sep=""),pch=21,cex=1.5,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1,cex.lab=1.2)
legend(x = "topright",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=0.7,pt.cex=1.5)
dev.off()


pdf("MDS_GATK_bias_DP.pdf",width=4.5,height = 4.5)

plot(x, y,
     xlab = "",
     ylab = "",yaxt="n",
     xaxt="n",main=main,col=color.inner,lwd=3,bg=paste(color.outer,"FF",sep=""),pch=21,cex=2,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.2),cex.axis=1.2,cex.lab=1.4, title(xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),line=1.2))
legend(x = "topright",legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 3,pch=21,cex=0.8,pt.cex=1.8)
axis(1, tck = -.01,mgp=c(3,0.2,0))
axis(2, tck = -.01,mgp=c(3,0.2,0))

dev.off()



ggsave(file="median_depth_samples.svg", plot=px, width=15, height=12)
dev.off()

pdf("median_depth_samples.pdf",width=5,height = 4, )
px
dev.off()


setwd("L:/MAVE/SCRATCH_BKP_2021/V7")
ggsave(file="corr_depthmax10_callers.svg", plot=p1, width=6, height=6)
dev.off()

pdf("corr_depthmax100_callers.pdf",width=5,height =5,useDingbats=FALSE)
p2
dev.off()



plot_grid(px,pmds,p1,p2,label_size = 28,labels = 'AUTO', hjust = 0, vjust = 1, scale= 0.8)

plot_grid(p1,p2,label_size = 28,labels = 'AUTO', hjust = 0, vjust = 1, scale= 0.8)
ppost<-p1+p2
setwd("L:/MAVE/SCRATCH_BKP_2021/V7")
ggsave(file="DP_phred_callers_poster.svg", plot=ppost, width=15, height=12)
dev.off()

pmds

###PLOT
FP<-px+pmds+p1+p2+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 28))

FP
#setwd("L:/MAVE/SCRATCH_BKP_2021/V7/")


pdf("DP_phred.pdf",width=15,height = 12)
FP
dev.off()


ggsave(file="DP_phred.svg", plot=FP, width=15, height=12)
dev.off()
