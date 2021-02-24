###geomorph nicoleti
                                               
library(geomorph)
library(StereoMorph)
library(RRPP)

####Analisis del cefalotorax#####
###################################

shapes.cefa <- readShapes("Shapes_cefa")
shapesGM.cefa <- readland.shapes(shapes.cefa)

#estimating missing
shapesGM.cefa <- estimate.missing(shapesGM.cefa, method = "TPS")
#no missing data

# Perform GPA Generalized Procrustes Analysis
Y.gpa.cefa <- gpagen(shapesGM.cefa, print.progress = FALSE)
plot(Y.gpa.cefa)
class(Y.gpa.cefa) #gpagen

coords <- Y.gpa.cefa$coords
class(coords)

csize <- Y.gpa.cefa$Csize
class(csize)
write.csv(csize, file = "csize_cefa.csv")

##### Data Pre-Processing

# Check for Outliers

plotOutliers(coords, inspect.outliers = T)

# Types of deformations

ref <- mshape(Y.gpa.cefa$coords)
par(mfrow=c(3,2))
plotRefToTarget(ref,Y.gpa.cefa$coords[,,27])
mtext("TPS")
plotRefToTarget(ref,Y.gpa.cefa$coords[,,27],mag=2.5)
mtext("TPS: 2.5X magnification")

plotRefToTarget(ref,Y.gpa.cefa$coords[,,27],method="vector",mag=3)
mtext("Vector Displacements")
plotRefToTarget(ref,Y.gpa.cefa$coords[,,27],gridPars=gridPar(pt.bg="red", link.col="green", pt.size = 1),
                method="vector",mag=3)
mtext("Vector Displacements: Other Options")

plotRefToTarget(ref,Y.gpa.cefa$coords[,,27],mag=2)  
mtext("Outline Deformation")
plotRefToTarget(ref,Y.gpa.cefa$coords[,,27],method="points")
mtext("OUtline Deformations Ref (gray) & and Tar (black)")
par(mfrow=c(1,1))

# Shape Predictions

# PCA-based
class(ref)
PCA <- plotTangentSpace(Y.gpa.cefa$coords)
PC <- PCA$pc.scores[,1]
preds <- shape.predictor(Y.gpa.cefa$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
par(mfrow=c(1,1))
plotRefToTarget(ref, preds$pred1)
mtext("PC1 - Min.")
plotRefToTarget(ref, preds$pred2)
mtext("PC1 - Max.")

# Regression-based
gdf <- geomorph.data.frame(Y.gpa.cefa)
nicoAllometry <- procD.lm(coords ~ log(Csize), data=gdf, print.progress = FALSE)
summary(nicoAllometry)
#Analysis of Variance, using Residual Randomization
#Permutation procedure: Randomization of null model residuals 
#Number of permutations: 1000 
#Estimation method: Ordinary Least Squares 
#Sums of Squares and Cross-products: Type I 
#Effect sizes (Z) based on F distributions
#           Df       SS       MS     Rsq      F      Z Pr(>F)   
#log(Csize)  1 0.0042601 0.0042601 0.17396 5.2647 3.1465  0.002 **
#Residuals  25 0.0202296 0.0008092 0.82604                        
#Total      26 0.0244897                                           
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Call: procD.lm(f1 = coords ~ log(Csize), data = gdf, print.progress = FALSE)

allom.plot <- plot(nicoAllometry, 
                   type = "regression", 
                   predictor = log(gdf$Csize),
                   reg.type ="PredLine") # make sure to have a predictor 
class(nicoAllometry)
preds <- shape.predictor(nicoAllometry$GM$fitted, x= allom.plot$RegScore, Intercept = FALSE, 
                         predmin = min(allom.plot$RegScore), 
                         predmax = max(allom.plot$RegScore)) 
class(preds)
plotRefToTarget(ref, preds$predmin, mag=3)
plotRefToTarget(ref, preds$predmax, mag=3)

# via picknplot.shape (more detail below)

picknplot.shape(allom.plot) 

# Group difference-based

data <- read.csv("data_geomorph_cepha.csv", header = TRUE, sep = ";")
View(data)

gdf2 <- geomorph.data.frame(ID = data$id, groups = data$Clado, k4 =data$K4, k5 = data$K5, COI = data$COI, BFD_gdi = data$BFD_gdi, size = data$csize, phy = tree)
attributes(gdf2)

# Single-Factor ANOVA

PCA1<- plotTangentSpace(Y.gpa.cefa$coords, groups = gdf2)

#clado
nicogm.anova <- procD.lm(coords ~ groups, data = gdf2, print.progress = FALSE)
anova(nicogm.anova)
#Analysis of Variance, using Residual Randomization
#Permutation procedure: Randomization of null model residuals 
#Number of permutations: 1000 
#Estimation method: Ordinary Least Squares 
#Sums of Squares and Cross-products: Type I 
#Effect sizes (Z) based on F distributions

#         Df        SS         MS     Rsq      F      Z Pr(>F)
#groups     2 0.0013913 0.00069566 0.05681 0.7228 -0.60639  0.724
#Residuals 24 0.0230984 0.00096243 0.94319                       
#Total     26 0.0244897                                        

#Call: procD.lm(f1 = coords ~ groups, data = gdf2, print.progress = FALSE)

#K=4
nicogm.anova1 <- procD.lm(coords ~ k4, data = gdf2, print.progress = FALSE)
anova(nicogm.anova1)
#          Df        SS        MS     Rsq      F      Z Pr(>F)
#k4         3 0.0030405 0.00101351 0.12416 1.0868 0.39116  0.339
#Residuals 23 0.0214491 0.00093257 0.87584                      
#Total     26 0.0244897   

#K=5
nicogm.anova2 <- procD.lm(coords ~ k5, data = gdf2, print.progress = FALSE)
anova(nicogm.anova2)
#          Df     SS        MS       Rsq      F      Z    Pr(>F)  
#k5         4 0.0042485 0.00106211 0.17348 1.1544 0.57319  0.293
#Residuals 22 0.0202412 0.00092006 0.82652                      
#Total     26 0.0244897     

#COI
nicogm.anova3 <- procD.lm(coords ~ COI, data = gdf2, print.progress = FALSE)
anova(nicogm.anova3)
#          Df        SS         MS     Rsq      F      Z Pr(>F)  
#COI         5 0.0058671 0.00117342 0.23957 1.3232 1.0643  0.156
#Residuals 21 0.0186226 0.00088679 0.76043                     
#Total     26 0.0244897  

#BFD_gdi
#nicogm.anova4 <- procD.lm(coords ~ BFD_gdi, data = gdf2, print.progress = FALSE)
#anova(nicogm.anova4)
#          Df     SS         MS     Rsq      F      Z     Pr(>F)  
#BFD_gdi    5 0.0070648 0.00141296 0.25011 1.4675 1.3901  0.092 .
#Residuals 22 0.0211820 0.00096282 0.74989                       
#Total     27 0.0282468  

##### 3: Principal Components Analysis (PCA)

plotTangentSpace(Y.gpa.cefa$coords, groups = data$Clado)

nic.raw <- gm.prcomp(Y.gpa.cefa$coords, groups = data$Clado)
nic.raw1 <- gm.prcomp(Y.gpa.cefa$coords, groups = data$COI)
gps <- as.factor(paste(data$id, data$región)) 
gps <- as.factor(data$Clado)
gps1 <- as.factor(data$K4) 
gps2 <- as.factor(data$K5)
gps3 <- as.factor(data$COI)
gps4 <- as.factor(data$X7_species)
plot(nic.raw)
plot(nic.raw1)
par(mar=c(2, 2, 2, 2))
par(mfrow=c(1,1))
plot(nic.raw, pch=21, cex = 1.5, bg = gps) 
plot(nic.raw, pch=21, cex = 1.5, bg = gps1) 
plot(nic.raw, pch=21, cex = 1.5, bg = gps2) 
plot(nic.raw, pch=21, cex = 1.5, bg = gps3)

#  Add things as desired using standard R plotting
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 37.68%", pos = 4, font = 2)
text(0, 0.95*par()$usr[4], labels = "PC2 - 17.13%", pos = 4, font = 2)
legend("bottomleft", pch=21, pt.bg = unique(gps3), legend = levels(gps3))


##############################################################################
##### 3: Phylogenetic Comparative Methods =======================================
# data input
library(geiger)
tree <- read.nexus('nico_cefa.tree')
plot(tree)
plot(ladderize(nicotree),edge.width=3)
axisPhylo(1)

### --------------------------------------------------------------------------------------------------

#PGLS regression
gdf3$
gdf3 <- geomorph.data.frame(Y.gpa.cefa, phy = tree) 
pgls.reg <- procD.pgls(coords ~ Csize, data=gdf3, phy=tree,iter = 999)
summary(pgls.reg)
#Analysis of Variance, using Residual Randomization
#Permutation procedure: Randomization of null model residuals 
#Number of permutations: 1000 
#Estimation method: Generalized Least-Squares (via OLS projection) 
#Sums of Squares and Cross-products: Type I 
#Effect sizes (Z) based on F distributions

#          Df      SS     MS       Rsq      F      Z   Pr(>F)   
#Csize         1 0.18437 0.184372 0.09613 2.6589 1.9829  0.029 *
#Residuals 25 1.73353 0.069341 0.90387                       
#Total     26 1.91791 
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Call: procD.lm(f1 = f1, iter = iter, seed = seed, RRPP = TRUE, SS.type = SS.type, effect.type = effect.type, int.first = int.first, Cov = Cov, data = data, print.progress = print.progress)


allom.plot <- plot(pgls.reg, type = "regression", predictor = gdf3$Csize, reg.type ="RegScore", pch=19, cex=1.5) # make sure to have a predictor 
#plots
preds <- shape.predictor(pgls.reg$GM$pgls.fitted, x= allom.plot$RegScore, Intercept = FALSE, 
                         predmin = min(allom.plot$RegScore), 
                         predmax = max(allom.plot$RegScore)) 
M <- mshape(coords)
plotRefToTarget(M, preds$predmin, mag=3, links = links)
plotRefToTarget(M, preds$predmax, mag=3, links = links)

### --------------------------------------------------------------------------------------------------
gdf4 <- geomorph.data.frame(Y.gpa.cefa, ID = data$id, clado = data$Clado, k4 =data$K4, k5 = data$K5, COI = data$COI, BFD_gdi = data$BFD_gdi, size = data$csize, phy = tree)
class(gdf4)

# PGLS ANOVA
#Phylogeentic ANOVA
#K=4
pgls.aov <- procD.pgls(coords~k4, effect.type = 'cohen', data=gdf4,
                       phy=tree, print.progress = FALSE)
summary(pgls.aov)
#Analysis of Variance, using Residual Randomization
#Permutation procedure: Randomization of null model residuals 
#Number of permutations: 1000 
#Estimation method: Generalized Least-Squares (via OLS projection) 
#Sums of Squares and Cross-products: Type I 
#Effect sizes (Z) based on F distributions

#          Df      SS      MS      Rsq      F       Z    Pr(>F)  
#k4          3 0.29531 0.098436 0.15397 1.3953 0.9183  0.179
#Residuals 23 1.62260 0.070548 0.84603                     
#Total     26 1.91791  

#clado
pgls.aov1 <- procD.pgls(coords~clado, effect.type = 'cohen', data=gdf4,
                       phy=tree, print.progress = FALSE)
summary(pgls.aov1)
#          Df       SS        MS     Rsq     F       Z Pr(>F)
#clado      2 0.17343 0.086713 0.09042 1.193 0.54834  0.288
#Residuals 24 1.74448 0.072687 0.90958                     
#Total     26 1.91791 

#K5
pgls.aov2 <- procD.pgls(f1 = coords ~ k5, effect.type = 'cohen', data=gdf4, phy=tree, print.progress = FALSE)
summary(pgls.aov2)
#          Df       SS        MS     Rsq      F       Z Pr(>F)
#k5         4 0.31729 0.079322 0.16543 1.0903 0.30825  0.377
#Residuals 22 1.60062 0.072755 0.83457                      
#Total     26 1.91791   

#COI
pgls.aov3 <- procD.pgls(f1 = coords ~ COI, effect.type = 'cohen', data=gdf4, phy=tree, print.progress = FALSE)
summary(pgls.aov3)
#          Df       SS         MS     Rsq      F       Z Pr(>F)
#COI        5 0.34576 0.069152 0.18028 0.9237 -0.16914  0.551
#Residuals 21 1.57215 0.074864 0.81972                       
#Total     26 1.91791  



#### Phylogenetic Signal
#The function estimates the degree of phylogenetic signal present in Procrustes shape variables for a given phylogeny. It is assumed that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) [e.g., with gpagen]. The degree of phylogenetic signal in data is estimated using the multivariate version of the K-statistic (Kmult: Adams 2014). This value evaluates the degree of phylogenetic signal in a dataset relative to what is expected under a Brownian motion model of evolution. For geometric morphometric data, the approach is a mathematical generalization of the Kappa statistic (Blomberg et al. 2003) appropriate for highly multivariate data (see Adams 2014).Significance testing is found by permuting the shape data among the tips of the phylogeny. In addition, a multivariate effect size describing the strength of the effect is estimated from the empirically-generated sampling distribution (see details in Adams and Collyer 2019). Values from these distributions are log-transformed prior to effect size estimation, to assure normally distributed data.
PS.shape <- physignal(coords, tree, print.progress = FALSE)
summary(PS.shape)
#Call:
#physignal(A = coords, phy = tree, print.progress = FALSE) 
#Observed Phylogenetic Signal (K): 0.123
#P-value: 0.04
#Effect Size: 1.6687
#Based on 1000 random permutations

plot(PS.shape) 


# Comparing Net Rates of Evolution
#ER<-compare.evol.rates(A=gdf4$coords, phy=nicotree,gp=gdf4$k4,iter=999,  method = 'permutation', print.progress = FALSE)
#summary(ER)
#plot(ER)  # COMPARISONS AMONG CLADES
#EMR <- compare.multi.evol.rates(A=gdf4$coords, phy=nicotree, gp=gps3, print.progress = FALSE)
#summary(EMR)
#plot(EMR)  # COMPARISONS AMONG TRAITS

##PhyloMorphoSpace
plotGMPhyloMorphoSpace(tree, coords)

##### 4: Allometry ===================================================================================

# Simple Allometry
fit <- procD.lm(coords ~ log(Csize), data=gdf4, iter=999, print.progress = FALSE)
anova(fit)
#Analysis of Variance, using Residual Randomization
#Permutation procedure: Randomization of null model residuals 
#Number of permutations: 1000 
#Estimation method: Ordinary Least Squares 
#Sums of Squares and Cross-products: Type I 
#Effect sizes (Z) based on F distributions

#           Df      SS       MS     Rsq      F      Z   Pr(>F)   
#log(Csize)  1 0.004120 0.004120 0.14586 4.4399 2.9925  0.001 **
#Residuals  26 0.024127 0.000928 0.85414                        
#Total      27 0.028247   

# Predline
plotAllometry(fit, size = gdf4$Csize, logsz = TRUE, method = "PredLine", pch = 19)

# RegScore
plotAllometry(fit, size = gdf4$Csize, logsz = TRUE, method = "RegScore", pch = 19)

# CAC
plotAllometry(fit, size = gdf4$Csize, logsz = TRUE, method = "CAC", pch = 19)

### --------------------------------------------------------------------------------------------------

# Group Allometry, including homogeneity of slopes test

#fit.unique <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=999, print.progress = FALSE)
#fit.common <- procD.lm(coords ~ Csize + species * site, data=gdf, iter=999, print.progress = FALSE)
#anova(fit.common, fit.unique)

# Because the unique slopes model was slightly better, it seems unwise to assume slopes
# are parallel and compare means.  However, the additional explained variation with unique slopes 
# was also quite small.
# Let's see what happens when we compare slopes.
# We can compare slopes with pairwise, just like means.
# Let's make sure the common slopes model is the null model.

### --------------------------------------------------------------------------------------------------



### --------------------------------------------------------------------------------------------------











#########################
#####QUELA###########
########################

##Digitizing Landmarks
#First, let´s see how one can digitize fixed landmarks. These correspond to well-defined anatomical points on your structure of interest. You will need a folder with the pictures you will digitize (Images, in our example) In this script, digitized data are stored in another folder, called Shapes

digitizeImages(image.file='quela_v', shapes.file='shapes',
              landmarks.ref=paste("LM", c(1:10), sep=""))

#To digitize for the first time:
#* Go to Landmarks tab to see list of landmarks
#* The first lm is in bold
#* In the image window, double-click on the desired anatomical location to digitize
#* The point is now selected, and can be edited, moved etc. (is highlighted in green)
#Select the next lm name and repeat
#Navigate between landmarks with “n” (next) and “p” (previous)
#Select lm and press “d” to delete
#Save before moving to the next picture (or select the automatic save option)
#NOTE: one has the ability to copy all landmarks to each image, and then move them to the desired locations. This can be useful so that landmarks remain in the correct order.
#Important detail! Any specimen with missing data (e.g., a broken specimen) may still be digitized. In this case, simply leave the missing landmarks BLANK. They will be treated as NA and can be estimated in subsequent data-processing steps.
#Scaling
#To scale landmark coordinates during superimposition and obtain accurate estimates of the relative size of the studied objects, you need to define the scaling of the images used, i.e. to indicate the correspondence of picture pixels with a “true” size measurement. For this:
#Go to the scaling tab
#Select the ruler points by double clicking (one can use >1 scaling intervals, and then an average is calculated)
#Introduce your reference interval and units
#Save

#readland.shapes. This is a geomorph function that reads in data from a class shapes object. (In geomorph, there are several readland.___ functions for different landmark file types.) Before using readland.shapes, a shapes object must be created (a Stereomorph function).
shapes <- readShapes("shapes")
shapesGM <- readland.shapes(shapes) # can change 10 to a different number

#estimating missing
shapesGM <- estimate.missing(shapesGM, method = "TPS")

# Perform GPA
Y.gpa <- gpagen(shapesGM, print.progress = FALSE)
plot(Y.gpa)
class(Y.gpa)

coords <- Y.gpa$coords
class(coords)

csize <- Y.gpa$Csize
class(csize)
write.csv(csize, file = "csize.csv")

#Later it might impress you that no arguments for gpagen are needed when reading data in with readland.shapes. Class geomorphShapes objects also have various attributes that can be exported for other uses.
class(shapesGM)
attributes(shapesGM)
shapesGM$fixed
shapesGM$curves
land <- shapesGM$landmarks
names(shapesGM$landmarks)
summary(shapesGM)

##### 3: Data Pre-Processing

# Check for Outliers

plotOutliers(coords, inspect.outliers = T)

# Fixed Angle
#jaw.fixed <- fixed.angle(Y.gpa$coords, art.pt=1, angle.pts.1 = 5, angle.pts.2 = 6, rot.pts = c(2,3,4,5))
#gpa.fixed <- gpagen(jaw.fixed, print.progress = FALSE)
#plotAllSpecimens(gpa.fixed$coords)

############################################
# Types of deformations

ref <- mshape(Y.gpa$coords)
par(mfrow=c(3,2))
plotRefToTarget(ref,Y.gpa$coords[,,47])
mtext("TPS")
plotRefToTarget(ref,Y.gpa$coords[,,47],mag=2.5)
mtext("TPS: 2.5X magnification")

plotRefToTarget(ref,Y.gpa$coords[,,47],method="vector",mag=3)
mtext("Vector Displacements")
plotRefToTarget(ref,Y.gpa$coords[,,47],gridPars=gridPar(pt.bg="red", link.col="green", pt.size = 1),
                method="vector",mag=3)
mtext("Vector Displacements: Other Options")

plotRefToTarget(ref,Y.gpa$coords[,,47],mag=2)  
mtext("Outline Deformation")
plotRefToTarget(ref,Y.gpa$coords[,,47],method="points")
mtext("OUtline Deformations Ref (gray) & and Tar (black)")
par(mfrow=c(1,1))

# Shape Predictions

# PCA-based
M <- mshape(Y.gpa$coords)
class(M)
PCA <- plotTangentSpace(Y.gpa$coords)
PC <- PCA$pc.scores[,1]
preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(M, preds$pred1)
mtext("PC1 - Min.")
plotRefToTarget(M, preds$pred2)
mtext("PC1 - Max.")

# Regression-based
gdf <- geomorph.data.frame(Y.gpa)
nicoAllometry <- procD.lm(coords ~ log(Csize), data=gdf, print.progress = FALSE)
summary(nicoAllometry)
#Analysis of Variance, using Residual Randomization
#Permutation procedure: Randomization of null model residuals 
#Number of permutations: 1000 
#Estimation method: Ordinary Least Squares 
#Sums of Squares and Cross-products: Type I 
#Effect sizes (Z) based on F distributions

#             Df       SS       MS     Rsq      F      Z  Pr(>F)   
#log(Csize)    1 0.004120 0.004120 0.14586 4.4399 2.9925  0.001 **
#  Residuals  26 0.024127 0.000928 0.85414                        
#Total        27 0.028247 

allom.plot <- plot(nicoAllometry, 
                   type = "regression", 
                   predictor = log(gdf$Csize),
                   reg.type ="PredLine") # make sure to have a predictor 
class(nicoAllometry)
preds <- shape.predictor(nicoAllometry$GM$fitted, x= allom.plot$RegScore, Intercept = FALSE, 
                         predmin = min(allom.plot$RegScore), 
                         predmax = max(allom.plot$RegScore)) 
class(preds)
plotRefToTarget(M, preds$predmin, mag=3)
plotRefToTarget(M, preds$predmax, mag=3)

# via picknplot.shape (more detail below)

picknplot.shape(allom.plot) 

# Group difference-based

data <- read.csv("data_geomorph.csv", header = TRUE, sep = ";")
gdf <- geomorph.data.frame(Y.gpa)
attributes(gdf)

gdf <- geomorph.data.frame(Y.gpa, indv = data$id, sp.cand = data$group, basin=data$basin, posi = data$geo, size = data$csize, phy = tree)
attributes(gdf)
# Single-Factor ANOVA

PCA <- plotTangentSpace(Y.gpa$coords, groups = gdf$species)

nicogm.anova <- procD.lm(coords ~ sp.cand, data = gdf, print.progress = FALSE)
anova(nicogm.anova)
#Analysis of Variance, using Residual Randomization
#Permutation procedure: Randomization of null model residuals 
#Number of permutations: 1000 
#Estimation method: Ordinary Least Squares 
#Sums of Squares and Cross-products: Type I 
#Effect sizes (Z) based on F distributions

#          Df        SS         MS    Rsq      F      Z  Pr(>F)  
#sp.cand    4 0.0063499 0.00158748 0.2248 1.6675 1.8087  0.036 *
#Residuals 23 0.0218969 0.00095204 0.7752                       
#Total     27 0.0282468                                         

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Call: procD.lm(f1 = coords ~ sp.cand, data = gdf, print.progress = FALSE)

X <- nicogm.anova$X

X # includes intercept; remove for better functioning 
X <- X[,-1]
symJord <- c(0,1,0) # design for P. Jordani in sympatry
alloJord <- c(0,0,0) # design for P. Jordani in allopatry
#preds <- shape.predictor(nicogm.anova$fitted, x = X, Intercept = TRUE, 
                        #symJord=symJord, alloJord=alloJord)

#plotRefToTarget(M, preds$symJord, links = plethodon$links, mag=2)
#plotRefToTarget(M, preds$alloJord, links = plethodon$links, mag=2)

# via picknplot.shape (more detail below)

plot.anova <- plot(nicogm.anova, type = "PC", pch = 21, 
                   bg = interaction(gdf$Csize, gdf$sp.cand), 
                   asp = 1)

plot.anova <- plot(nicogm.anova, type = "PC", pch = 21, 
                   bg = interaction(gdf$Csize, gdf$posi), 
                   asp = 1)



##### 3: Principal Components Analysis (PCA)

plotTangentSpace(Y.gpa$coords, groups = data$group)

nic.raw <- gm.prcomp(Y.gpa$coords)

gps <- as.factor(paste(data$id, data$group)) 
gps <- as.factor(data$group)
gps <- as.factor(data$geo) 
plot(nic.raw)
par(mar=c(2, 2, 2, 2))
plot(nic.raw, pch=22, cex = 1.5, bg = gps) 
#  Add things as desired using standard R plotting
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", pos = 4, font = 2)
text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2)
legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps))


##### 3: Phylogenetic Comparative Methods ============================================================

# data input
library(geiger)
tree <- read.nexus('nico_cefa.tree')
plot(tree)
class(tree)

#dat <- read.csv('Data/svl.csv', header=TRUE, row.names=1)
#svl <-dat[,1]; names(svl) <- rownames(dat)
#shape <- readland.tps('Data/headshape.tps',specID = "ID",warnmsg = FALSE)

#match.data <- treedata(plethtree,svl)  

#plethgps <- read.csv('Data/Gps.csv',header=TRUE, row.names=1)
#plethgps <- plethgps[match(dimnames(shape)[[3]],rownames(plethgps)),]
#elev <- as.factor(plethgps$ElevGp); names(elev) <- rownames(plethgps)

#gdf <- geomorph.data.frame(shape=shape, svl=svl,elev = elev, plethtree=plethtree)
#links <- matrix(c(4,3,2,1,1,6,7,8,9,10,1,1,11,5,5,4,2,3,7,8,9,10,11,9,10,1),
#                ncol=2,byrow=FALSE)
plot(ladderize(tree),edge.width=3)
axisPhylo(1)

### --------------------------------------------------------------------------------------------------

# PGLS regression
pgls.reg <- procD.pgls(f1 = coords ~ Csize, effect.type = 'cohen', data=gdf, phy=tree, print.progress = FALSE)
summary(pgls.reg)
allom.plot <- plot(pgls.reg, type = "regression", predictor = gdf$sais, reg.type ="RegScore", pch=19, cex=1.5) # make sure to have a predictor 
#plots
preds <- shape.predictor(pgls.reg$GM$pgls.fitted, x= allom.plot$RegScore, Intercept = FALSE, 
                         predmin = min(allom.plot$RegScore), 
                         predmax = max(allom.plot$RegScore)) 
M <- mshape(coords)
plotRefToTarget(M, preds$predmin, mag=3, links = links)
plotRefToTarget(M, preds$predmax, mag=3, links = links)

### --------------------------------------------------------------------------------------------------

# PGLS ANOVA
pgls.aov <- procD.pgls(f1 = coords~posi, effect.type = 'cohen', data=gdf,
                       phy=nicotree, print.progress = FALSE)
summary(pgls.aov)

pc.plot <- plot(pgls.aov,type = "PC", pch=21, cex=1.5,bg=gdf$posi)

shapeHulls(pc.plot, groups = gdf$posi, 
           group.cols = c("red", "black"),
           group.lwd = rep(1, 2), group.lty = c(2, 1))
legend("topright", levels(gdf$posi), 
       col = c("black", "red"),
       lwd = rep(1,2), lty = c(2, 1))

# plots
Low <- c(1) # design for low elevation
High <- c(0) # design for high elevation
preds <- shape.predictor(arrayspecs(pgls.aov$pgls.fitted, 0,2), x= pgls.aov$X[,-1],
                         Intercept = TRUE, Low = Low, High = High)   
par(mfrow=c(1,2)) 
plotRefToTarget(M, preds$Low, mag=2, links=links)
mtext("Low Elevation")
plotRefToTarget(M, preds$High, mag=2, links=links)
mtext("High Elevation")
par(mfrow=c(1,1)) 

# via picknplot.shape
picknplot.shape(pc.plot)

### --------------------------------------------------------------------------------------------------
nicotree$edge.length
# Phylogenetic PLS
land.gps<-c("A","A","A","A","A","B","B","B","B","B","B")
PLS.Y <- phylo.integration(A = gdf$sais, partition.gp = gps, 
                           phy= nicotree, print.progress = FALSE)
summary(PLS.Y)
pls.plot <- plot(PLS.Y)
picknplot.shape(pls.plot, mag = 2)

# Phylogenetic Ordination
plot.res <- gm.prcomp(coords,phy=tree, data=gdf4)
ord.plot <- plot(plot.res,phylo = TRUE, pch=21, bg="red", cex=1.5)
picknplot.shape(ord.plot, mag = 2)

# Phylogenetic Signal
PS.shape <- physignal(coords, nicotree, print.progress = FALSE)
summary(PS.shape)
plot(PS.shape)

# Comparing Net Rates of Evolution
ER<-compare.evol.rates(A=gdf4$coords, phy=tree,gp=gdf4$k5,iter=999, 
                       method = 'permutation', print.progress = FALSE)
summary(ER)
plot(ER)  # COMPARISONS AMONG CLADES

EMR <- compare.multi.evol.rates(A=gdf$coords, phy=tree, gp=gps1, print.progress = FALSE)
summary(EMR)
plot(EMR)  # COMPARISONS AMONG TRAITS

##### 4: Allometry ===================================================================================

# Simple Allometry
data(plethodon) 
Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment  
gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
                           species = plethodon$species) 
fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=999, print.progress = FALSE)
anova(fit)

# Predline
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", pch = 19)

# RegScore
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "RegScore", pch = 19)

# CAC
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "CAC", pch = 19)

### --------------------------------------------------------------------------------------------------

# Group Allometry, including homogeneity of slopes test

fit.unique <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=999, print.progress = FALSE)
fit.common <- procD.lm(coords ~ Csize + species * site, data=gdf, iter=999, print.progress = FALSE)
anova(fit.common, fit.unique)

# Because the unique slopes model was slightly better, it seems unwise to assume slopes
# are parallel and compare means.  However, the additional explained variation with unique slopes 
# was also quite small.
# Let's see what happens when we compare slopes.
# We can compare slopes with pairwise, just like means.
# Let's make sure the common slopes model is the null model.

### --------------------------------------------------------------------------------------------------

# Pairwise comparisons
slope.pw <- pairwise(fit.unique, fit.null = fit.common, 
                     groups = interaction(gdf$species, gdf$site),
                     covariate = gdf$Csize)
summary(slope.pw, test.type = "VC", angle.type = "deg") # angular differences
summary(slope.pw, test.type = "dist", angle.type = "deg") # amount of shape change differences

# Conclusion: some slight differences in angles between slopes
# Note that the UCL angles are quite large - usually an indication of 
# small size ranges.

### --------------------------------------------------------------------------------------------------
# Plots

# Predline
plotAllometry(fit.unique, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
              pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))

# RegScore
plotAllometry(fit.unique, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
              pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))


# Size-Shape Space
pc.plot <- plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "size.shape", 
                         pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
summary(pc.plot$size.shape.PCA)


##### 4: PickNPlot Shapes in Real Time (more detail here)

nico.pca <- gm.prcomp(Y.gpa$coords)

pca.plot <- plot(nico.pca)
picknplot.shape(pca.plot) 

picknplot.shape(plot(nico.pca), method = "points", mag = 3)

##### 5: 3D Warping

scallops <- readland.tps("Data/scallops for viz.tps", specID = "ID")
ref <- mshape(scallops)
refmesh <- warpRefMesh(read.ply("Data/glyp02L.ply"), 
                       scallops[,,1], ref, color=NULL, centered=T)
plotTangentSpace(scallops, axis1 = 1, axis2 = 2, warpgrids=T, mesh= refmesh)

##### 6: Two-Block Partial Least Squares (PLS)

data(pupfish)
Y.gpa <- gpagen(pupfish$coords, print.progress = FALSE)
plotAllSpecimens(Y.gpa$coords)
shape <- Y.gpa$coords
headland <- c(4, 10:17, 39:56)

PLS <- two.b.pls(shape[headland,,],shape[-headland,,], iter=999, print.progress = FALSE)
summary(PLS)
pls.plot <- plot(PLS)

## PLS shape predictions
preds <- shape.predictor(shape[headland,,], two.d.array(shape[-headland,,]), Intercept = FALSE,
                         method = "PLS", pred1 = -0.2, pred2 = 0.2) # using PLS plot as a guide

M <- mshape(shape[headland,,])
plotRefToTarget(M, -1*preds$pred1, mag=3)
plotRefToTarget(M, preds$pred2, mag=3)

# via picknplot.shape (more detail above)

picknplot.shape(pls.plot, mag = 3) 

##### 7: Regression

pupfish$logSize <- log(pupfish$CS)  #add logCS to geomorph data frame
fit <- procD.lm(coords ~ logSize, data = pupfish, print.progress = FALSE)
anova(fit)

plot(fit)
plot(fit, type = "regression", reg.type = "PredLine", predictor = pupfish$logSize, pch=21, bg="red")
plot(fit, type = "regression", reg.type = "RegScore", predictor = pupfish$logSize, pch=21, bg="red")

## Regression predictions
allom.plot <- plot(fit, 
                   type = "regression", 
                   predictor = pupfish$logSize,
                   reg.type ="RegScore") # make sure to have a predictor 

preds <- shape.predictor(allom.plot$GM$fitted, x= allom.plot$RegScore, Intercept = FALSE, 
                         predmin = min(allom.plot$RegScore), 
                         predmax = max(allom.plot$RegScore)) 
M <- mshape(pupfish$coords)
plotRefToTarget(M, preds$predmin, mag=3)
plotRefToTarget(M, preds$predmax, mag=3)


# via picknplot.shape (more detail above)

picknplot.shape(allom.plot, mag = 3) 

