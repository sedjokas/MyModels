# COde R pour recuperer les plots des data recu de GAMA
library(readr)
library(dplyr)

#For S=1000
#SIR si on regroupe quelques compartiments
datasir<-read_csv(file = "~/MyGAMAPhD/sante4Dev/TBC/results/datasir_.csv")
g_range <- range(0, datasir[,1])
k<- range(0,1100)
plot(datasir[,c(1,2)], type="l", col="green", axes=T, ann=T, xlab="Time (Years) ", ylab="Number of S(t), I(t) and R(t)   x100", ylim=k, , xlim= g_range)
lines(datasir[,c(1,3)], type="l", col="red")
lines(datasir[,c(1,4)], type="l", col="blue")

box()
legend(290, 1100, c("S(t)", "I(t)", "R(t)"), cex=0.8, 
   col=c("green", "red", "blue"), lty=1
);

# ALL si on garde l'equation telle
datasir<-read_csv(file = "~/MyGAMAPhD/sante4Dev/TBC/results/dataall.csv")
g_range <- range(0, datasir[,1])
k<- range(0,1100)
plot(datasir[,c(1,2)], type="l", col="green", axes=T, ann=T, xlab="Time (Years) ", ylab="Number of S(t), I(t), Le(t) , Lf(t) , T(t) , K(t) , R1(t) and R2(t)  x100", ylim=k, , xlim= g_range)
lines(datasir[,c(1,3)], type="l", col="red")
lines(datasir[,c(1,4)], type="l", col="orange")
lines(datasir[,c(1,5)], type="l", col="yellow")
lines(datasir[,c(1,6)], type="l", col="gray")
lines(datasir[,c(1,7)], type="l", col="magenta")
lines(datasir[,c(1,8)], type="l", col="blue")
lines(datasir[,c(1,9)], type="l", col="#77B5FE")

box()
legend(290, 1100, c("S(t)", "I(t)", "Le(t)" , "Lf(t)" , "T(t)" , "K(t)" , "R1(t)", "R2(t)"), cex=0.8, 
   col=c("green", "red", "orange", "yellow", "gray", "magenta", "blue", "#77B5FE"), lty=1);
