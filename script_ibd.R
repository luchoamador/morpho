####Isolation by distance analysis

library("adegenet")
library("ade4")
library("MASS")

####Isolation by distance analysis

#1
my.genind <-  read.structure(file="myfile.str")
allele.counts <- tab(my.genind, freq=TRUE, NA.method="mean")

#Distance Matrix Computation
Dgen <- dist(allele.counts)

#2
my.coords <- read.csv(file = "coord.csv", header = TRUE)
Dgeo <- dist(my.coords)

##Isolation by distance
ibd <- mantel.randtest(matriz, Dgeo)

plot(ibd)
#there is or is not a clear isolation by distance pattern.

###Cline or distant patches?
#The correlation between genetic and geographic distances can occur under a range of different biological scenarios. Classical IBD would result in continuous clines of genetic differentiation and cause such correlation. However, distant and differentiated populations would also result in such a pattern. These are slightly different processes and we would like to be able to disentangle them. A very simple first approach is simply plotting both distances:

plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo), col="red",lty=2)

#Most of the time, simple scatterplots fail to provide a good picture of the data as the density of points in the scatterplot is badly displayed. Colors can be used to provide better (and prettier) plots. Local density is measured using a 2-dimensional kernel density estimation (kde2d), and the results are displayed using image; colorRampPalette is used to generate a customized color palette:


dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")