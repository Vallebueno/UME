#Paper UME  Figure5
Sys.setenv(LANG = "en")

library("data.table")
library("ggplot2")
library("ggpmisc")
library("patchwork")
#install.packages("ggridges")
library("ggridges")
library("colourvalues")

###names(fit$points[which(fit$points[,2]>.1),]) ## find the weird samples




setwd("L:/MAVE/Experiments/DB_COF_alldata")






file='UME_production_2023_08_24.vcf.gz.smp.vcf_IBS.dist.txt'

### SFS
par(mfrow = c(1, 1), mar=c(3.4,3.3,0.7,0.5), mgp = c(1.7, 0.5, 0))

#par(mfrow = c(1, 1), mar=c(3,3.5,0.5,0.6))

colrs<-colour_values(1:17, palette="matlab_like", alpha=255)
#colrs<-colour_values(1:9, palette="greys", alpha=255)

file='UME_production_2023_08_24.vcf.gz.smp.vcf.SFS'
dHM <- na.omit(fread(file))

#title(alg, line = 0.6, cex.main = 2)
#par(mfrow = c(5, 1), mar=c(2.7,3.5,0.5,0.6), mgp = c(1.7, 0.5, 0))

#plot(density(dHM$SFS),col=colrs[1],lwd=2,xlab="MAF",fo ylab="Density", main="", cex.axis=1.2, cex.lab=1.2,ylim=c(0,15),cex.main=1.2)
line<-par(lwd=3)
hist(dHM$SFS,freq = FALSE,col=colrs[1],lwd=2,xlab="MAF", ylab="Density", main="", cex.axis=1.2, cex.lab=1.2,ylim=c(0,13.5),cex.main=1.2)


smallest_number <- min(dHM$SFS)
threshold <-0.0005055613
count <- length(dHM$SFS[dHM$SFS < threshold])
count
list_size <- length(dHM$SFS)
aver <- median(dHM$SFS)

(count/list_size)*100

SFS<-recordPlot()

#Metrics distribution
#den=density(dHM$SFS)
#den$x[which(den$y==max(den$y))]
#max(den$y)


####SITES summary

file='UME_production_2023_08_24.vcf.gz.smp.vcf_SiteSummary.txt'


data <- fread(file)
min(data$`Proportion Missing`) #describe in text

mean(data$`Proportion Missing`)


hist(data$`Proportion Missing`,freq=FALSE,col=colrs[1],lwd=2,xlab="Proportion missing per site", ylab="Density", main="", cex.axis=1.2, cex.lab=1.2,ylim=c(0,8),cex.main=1.2)
PMSNP<-recordPlot()

min(data$`Proportion Heterozygous`) #describe in text

mean(data$`Proportion Heterozygous`)


hist(na.omit(data$`Proportion Heterozygous`),freq=FALSE,col=colrs[1],lwd=2,xlab="Heterozygosity per site", ylab="Density", main="", cex.axis=1.2, cex.lab=1.2,ylim=c(0,7),cex.main=1.2)

HETSNP<-recordPlot()



### Taxa summary

file='UME_production_2023_08_24.vcf.gz.smp.vcf_TaxaSummary.txt'

datax <- fread(file)

hist(na.omit(datax$`Proportion Heterozygous`),freq=FALSE,col=colrs[1],lwd=2,xlab="Heterozygosity per taxon", ylab="Density", main="", cex.axis=1.2, cex.lab=1.2,ylim=c(0,5),cex.main=1.2)

HETAX<-recordPlot()


max(datax$`Proportion Missing`)
mean(datax$`Proportion Missing`)

hist(datax$`Proportion Missing`,freq=FALSE,col=colrs[1],lwd=2,xlab="Proportion missing per taxon", ylab="Density", main="", cex.axis=1.2, cex.lab=1.2,ylim=c(0,4),cex.main=1.2)
PMTAX<-recordPlot()


threshold <-0.0005055613
count <- length(dHM$SFS[datax$`Proportion Missing`])







##########MDS 1
par(mfrow = c(1, 1), mar=c(3.2,4.2,2,1),bg=NA, mgp = c(2, 0.5, 0))

#par(mfrow = c(1, 1), mar=c(2.7,3.5,0.5,0.6), mgp = c(1.7, 0.5, 0))

file='UME_production_2023_08_24.vcf.gz.smp.vcf_IBS.dist.txt'

colrf='Colors_ALL_August.txt'

cuan=1
distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
colnames(distancesAll) <- rownames(distancesAll)
col <- read.table(colrf,header=T,sep="\t",as.is=T, fill=TRUE)
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
outerClass  <- "group2"
innerClass <- "group2"
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.outer.legend <- c(col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]
i <-1
j <- i+1
x <- fit$points[,i]
y <- fit$points[,j]

color.outer.legend <- c(col[match(legend,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])


plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.outer,lwd=2,bg=color.inner,pch=21,cex=0.7,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*0.02),cex.axis=1.2, cex.lab=1.2,cex.main=0.1)
#legend("center",legend = legend,col=color.outer.legend,pt.lwd = 0,pch=16,cex=0.8,pt.cex=1,bty="n",y.intersp=0.9)
#legend(-0.17,-0.01,legend = legend,col=color.inner.legend,pt.bg=color.outer.legend,pt.lwd = 2,pch=21,cex=0.75,pt.cex=1,bty="n",y.intersp=2)
pmds1<-recordPlot()

which(fit$points[,2]>.1)


par(mfrow = c(1, 1), mar=c(0.01,1,0.01,0.01),bg=NA,oma = c(0, 0, 0, 0))
plot.new()
legend <- c("Nigeria", "Thailand", "South Africa","China", "Canada", "USA", "Mexico", "Guatemala", "Colombia", "Venezuela", "Ecuador", "Bolivia", "Peru", "Paraguay", "Brazil", "Chile", "Zea mays ssp. parviglumis", "Zea mays ssp. mexicana", "Zea mays ssp. huehuet.", "Zea luxurians", "Zea diploperennis", "Zea perennis", "Zea nicaraguensis", "Tripsacum dactyloides")
color.outer.legend <- c(col[match(legend,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend("left",legend = legend,col=color.outer.legend, ncol=3,pt.lwd = 0,pch=16,cex=1,pt.cex=1.2,bty="n",y.intersp=2.5,x.intersp=1)

pleg<-recordPlot()
pleg



##########MDS 2
par(mfrow = c(1, 1), mar=c(3.2,4.2,2,1),bg=NA, mgp = c(2, 0.5, 0))
file='UME_production_2023_08_24.vcf.gz.smp.vcf_IBS.dist.txt'

colrf='Colors_Nteo_August.txt'

cuan=1
distancesAll <- read.table(file,sep="\t",header=F,skip=5,row.names=1,as.is=T)
colnames(distancesAll) <- rownames(distancesAll)
col <- read.table(colrf,header=T,sep="\t",as.is=T, fill=TRUE)
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
outerClass  <- "group2"
innerClass <- "group2"
legend.outer <- levels(factor(col[match(rownames(fit$points[]),col$Taxon,),outerClass],levels = unique(col[,outerClass])))
color.outer.legend <- c(col[match(legend.outer,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
legend <- c(legend.outer)
color.inner <- col[match(rownames(fit$points[]),col$Taxon,),paste(innerClass,"color",sep="_")]
color.outer <- col[match(rownames(fit$points[]),col$Taxon,),paste(outerClass,"color",sep="_")]
i <-1
j <- i+1
x <- fit$points[,i]
y <- fit$points[,j]
legend <- c("Nigeria", "Thailand", "South Africa","China", "Canada", "USA", "Mexico", "Guatemala", "Colombia", "Venezuela", "Ecuador", "Bolivia", "Peru", "Paraguay", "Brazil", "Chile")
color.outer.legend <- c(col[match(legend,as.vector(col[,outerClass])),paste(outerClass,"color",sep="_")])
plot(x, y, xlab=paste("Coordinate ",i," (",round((fit$eig[i]/sum(fit$eig))*100,1),"%)",sep=""), ylab=paste("Coordinate ",j," (",round((fit$eig[j]/sum(fit$eig))*100,1),"%)",sep=""),main="",col=color.outer,lwd=2,bg=color.inner,pch=21,cex=0.5,xlim=c(min(fit$points[,1]),max(fit$points[,1])+(max(fit$points[,1])-min(fit$points[,1]))*.01),cex.axis=1.2, cex.lab=1.2,cex.main=0.1)
#legend(0.04,0.08,legend = legend,col=color.outer.legend,pt.lwd = 0,pch=16,cex=0.9,pt.cex=1,bty="n",y.intersp=1.1)
pmds2<-recordPlot()






### COW ploting 
#install.packages("gridBase")
library("gridBase")
library("cowplot")

FIGA<-plot_grid(pmds1,pleg,pmds2,
          nrow = 3,
          labels =  c('A','','B'),
          label_size = 33,
          vjust = 1.3,
          hjust = 0,
          align = "v"
)
#FIGA

FIGB<-plot_grid(SFS, PMSNP, HETSNP, HETAX, PMTAX,
                nrow = 5,
                ncol = 1,
                labels = c('C','D','E','F','G'),
                label_size = 33,
                vjust = 1.3,
                hjust = 0.5,
                align = "v"
)
#FIGB


FIGF<-plot_grid(FIGA,FIGB,label_size = 33,labels = "")
#FIGF


#marginstext


#pdf("Vallebueno2022_Fig_5.pdf",width=6.69,height = 6.5)
pdf("Vallebueno2023_1_Fig_5.pdf",width=10,height = 10)
FIGF
dev.off()
