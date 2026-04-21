library("colourvalues")
library("data.table")
setwd("L:/MAVE/TEST")

par(mfrow = c(2, 1), mar=c(4,5,2,1))

###### Not correlated
data <- data.frame(x=rnorm(100,20,0.7), y=rnorm(100,30,2))


plot(x=data$x,y=data$y,  xlab="Annual temperature °C",ylab="ring width mm", main="Not correlated")
abline(lm(data$y~data$x), col="red", lwd=4.0)


data <- na.omit(fread("Ancient_samples_info.txt"))


###### correlated

data$rem<-data$PC1*100+24

data$rrex<-(data$mean_genome_coverage/2.7)
plot(x=data$rem,y=data$rrex, xlab="Annual temperature °C",ylab="ring width mm", main="Correlated")
abline(lm(data$rrex~data$rem), col="red", lwd=4.0)





##########Genetic Distance



par(mfrow = c(2, 1), mar=c(4,5,2,1))

###### Not correlated
data <- data.frame(x=rnorm(100,0.5,0.1), y=rnorm(100,30,2))


plot(x=data$x,y=data$y,  xlab="Genetic Distance",ylab="Average ring width mm", main="Not correlated")
abline(lm(data$y~data$x), col="red", lwd=4.0)


data <- na.omit(fread("Ancient_samples_info.txt"))


###### correlated

data$rem<-(data$PC1*9)+1

data$rrex<-(data$mean_genome_depth/2.7)
plot(x=data$rem,y=data$rrex, xlab="Genetic Distance",ylab="Average ring width mm", main="Correlated")
abline(lm(data$rrex~data$rem), col="red", lwd=4.0)




#########picea fake data
par(mfrow = c(2, 1), mar=c(4,5,2,1))
data <- na.omit(fread("Picea_fake_data.txt"))
plot(x=data$Temp,y=data$mm, xlab="Temperature in C",ylab="ring width mm", main="Correlated")
abline(lm(data$mm~data$Temp), col="red", lwd=4.0)



new<-rnorm(12,20,0.7)

data$Temp2<- c(new)


plot(x=data$Temp2,y=data$mm, xlab="Temperature in C",ylab="ring width mm", main="Not Correlated")
abline(lm(data$mm~data$Temp2), col="red", lwd=4.0)



plot(x=data$Year,y=data$Temp, xlab="Year",ylab="Temperature in C", main="Correlated")
plot(x=data$Year,y=data$Temp2, xlab="Year",ylab="Temperature in C", main="Not Correlated")
