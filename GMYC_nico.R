###Species delimitation with GMYC
###Luis Amador

install.packages("splits", repos="http://R-Forge.R-project.org")

library(splits)
library(rncl)

###Primer Ejercicio con seqs incluido el outgroup y priors en BEAST como yule calibration y reloj strict
tre <- read.tree(file = "tree_calibrated.newick")
str(tre)
class(tre)
summary(tre)
tre
is.ultrametric(tre)
result <- gmyc(tre)
result

summary(result) 
#Result of GMYC species delimitation
#method:	single
#likelihood of null model:	477.0851
#maximum likelihood of GMYC model:	487.9154
#likelihood ratio:	21.66061
#result of LR test:	1.979058e-05***
#number of ML clusters:	8
#confidence interval:	6-16
#number of ML entities:	19
#confidence interval:	17-29
#threshold time:	-0.05833242

plot(result)
###los resultados son muy splitters sobredelimitan las especies... en BEAST estimation use birth death model, reloj stricto. será mejor usar los haplotipos para el analisis o incluir el outgroup.

##Segundo ejercicio con los mismso datos y parametros pero con umbral multiple de GMYC
result2 <- gmyc(tre, method = "multiple")
summary(result2)
plot(result2)
###Igualmente existe sobre estimación en el número de especies putativas.


################################################################


####bGMYC####
install.packages("bGMYC", repos = "https://pedrosenna.github.io/drat/")

library(bGMYC)

tre <- force.ultrametric(tre)

bgmyc.singlephy(tre, mcmc = 100000, burnin = 1, thinning = 10, t1=2, t2=88, start = c(1,1,25)) -> result.single

bgmyc.spec(result.single)->result.espec

spec.probmat(result.single)->result.probmat

str(result.probmat)

result.probmat

checkrates(result.single)->bgmyctasas

plot(result.probmat, tre)

plot(bgmyctasas)

bgmyctasas


library(phytools)
require(maps)

obj<-ltt(tre,plot=FALSE)

plot(obj,log="y",log.lineages=FALSE,main="lineage through time plot")
tree<-pbtree(b=1,d=0.25,t=4)
obj<-ltt(tree,gamma=FALSE)
obj

library(paleotree)

ltt.plot(tre, xlab = "Time", ylab = "N",
         backward = TRUE, tol = 1e-6)
ltt.lines(tre, backward = TRUE, tol = 1e-6, lty=2)

mltt.plot(tre, dcol = TRUE, dlty = T, legend = F, xlab = "Time", ylab = "N")

bd.cray <- birthdeath(tre)
bd.cray


####otro ejemplo de GMYC
library(ape)
library(splits)

Name.tr<-read.tree("*.nwk", tree.names = NULL)
##read fully-resolved and ultrametric input tree constructed with haplotype data 

test1 <- gmyc(Name.tr, method="single", interval=c(0, 10))
##execute single threshold analysis

summary(test1)  
##show summary results

plot(test1)
##plot results

test2 <- gmyc(Name.tr, method="multiple", interval= c(0, 10))
##execute multiple threshold analysis

summary(test2)
##show summary results

plot(test2)
##plot results

compare(test1,test2) 
##compare multiple vs single threshold models




