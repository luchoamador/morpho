###morpho analysis P. nicoleti###

library(ggpubr)
library(moments)

nico <- read.table("data_morpho.csv", header = TRUE, sep = ";")
nico
View(nico)
str(nico)

##Excluding (dropping) variables
myvars <- names(nico) %in% c("Intersex", "Other.aspects")
nico_morpholin <- nico[!myvars]
str(nico_morpholin)
View(nico_morpholin)

#Skewness is a measure of symmetry for a distribution. The value can be positive, negative or undefined.
#The skewness coefficient can be computed using the moments R packages:
skewness(nico_morpholin$TL, na.rm = TRUE)
#[1] 0.2413087
skewness(nico_morpholin$TL + nico_morpholin$CL + nico_morpholin$Cel, na.rm = TRUE)
#[1] 0.3094541
skewness(nico_morpholin$RL + nico_morpholin$RW, na.rm = TRUE)
#[1] 1.000835
skewness(nico_morpholin$DL + nico_morpholin$PrW + nico_morpholin$PrL, na.rm = TRUE)
#[1] 0.2638375
skewness(nico_morpholin$x + nico_morpholin$y + nico_morpholin$z, na.rm = TRUE)
#[1] 0.04040531
skewness(nico_morpholin$AreL + nico_morpholin$AreW, na.rm = TRUE)
#[1] -0.01302294
skewness(nico_morpholin$AreL, na.rm = TRUE)
#[1] 0.05646804
skewness(nico_morpholin$AreW, na.rm = TRUE)
#[1] 0.3513903


###Kernel density plot of total lenght of samples, kernel density estimation is a nonparametric method for estimating the probability density function of a random variable. Kernel density
###plots can be an effective way to view the distribution of a continuous variable.

#todos
all <- density(nico_morpholin$TL + nico_morpholin$CL + nico_morpholin$Cel + nico_morpholin$RL + nico_morpholin$RW + nico_morpholin$PrL + nico_morpholin$DL + nico_morpholin$PrW + nico_morpholin$x + nico_morpholin$y + nico_morpholin$z + nico_morpholin$AreW + nico_morpholin$AreL + nico_morpholin$AW)
plot(all, main = "Kernel Density of all characters")
polygon(all, col = "red", border = "black")
rug(nico_morpholin$TL + nico_morpholin$CL + nico_morpholin$Cel + nico_morpholin$RL + nico_morpholin$RW + nico_morpholin$PrL + nico_morpholin$DL + nico_morpholin$PrW + nico_morpholin$x + nico_morpholin$y + nico_morpholin$z + nico_morpholin$AreW + nico_morpholin$AreL + nico_morpholin$AW, col = "brown")

#Largo total
dl <- density(nico_morpholin$TL + nico_morpholin$CL + nico_morpholin$Cel)
plot(dl, main = "Kernel Density of Total Lenght")
polygon(dl, col = "red", border = "black")
rug(nico_morpholin$TL + nico_morpholin$CL + nico_morpholin$Cel, col = "brown")

dl.1 <- density(nico_morpholin$TL)
plot(dl.1, main = "Kernel Density of TL")
polygon(dl.1, col = "red", border = "black")
rug(nico_morpholin$TL, col = "brown")

dl.2 <- density(nico_morpholin$CL)
plot(dl.2, main = "Kernel Density of CL")
polygon(dl.2, col = "red", border = "black")
rug(nico_morpholin$CL, col = "brown")

dl.3 <- density(nico_morpholin$Cel)
plot(dl.3, main = "Kernel Density of Cel")
polygon(dl.3, col = "red", border = "black")
rug(nico_morpholin$Cel, col = "brown")

#Rostro
er <- density(nico_morpholin$RL + nico_morpholin$RW)
plot(er, main = "Kernel Density of rostrum")
polygon(er, col = "red", border = "black")
rug(nico_morpholin$RL + nico_morpholin$RW, col = "brown")

er.1 <- density(nico_morpholin$RL)
plot(er.1, main = "Kernel Density of RL")
polygon(er.1, col = "red", border = "black")
rug(nico_morpholin$RL, col = "brown")

er.2 <- density(nico_morpholin$RW)
plot(er.2, main = "Kernel Density of RW")
polygon(er.2, col = "red", border = "black")
rug(nico_morpholin$RW, col = "brown")

#quela
nq <- density(nico_morpholin$PrL + nico_morpholin$DL + nico_morpholin$PrW)
plot(nq, main = "Kernel Density of chelae")
polygon(nq, col = "red", border = "black")
rug(nico_morpholin$PrL + nico_morpholin$DL + nico_morpholin$PrW, col = "brown")

nq.1 <- density(nico_morpholin$PrL)
plot(nq.1, main = "Kernel Density of PrL")
polygon(nq.1, col = "red", border = "black")
rug(nico_morpholin$PrL, col = "brown")

nq.2 <- density(nico_morpholin$DL)
plot(nq.2, main = "Kernel Density of DL")
polygon(nq.2, col = "red", border = "black")
rug(nico_morpholin$DL, col = "brown")

nq.3 <- density(nico_morpholin$PrW)
plot(nq.3, main = "Kernel Density of PrW")
polygon(nq.3, col = "red", border = "black")
rug(nico_morpholin$PrW, col = "brown")

#Second pleura
sp <- density(nico_morpholin$x + nico_morpholin$y + nico_morpholin$z)
plot(sp, main = "Kernel Density of Second pleura")
polygon(sp, col = "red", border = "black")
rug(nico_morpholin$x + nico_morpholin$y + nico_morpholin$z, col = "brown")

sp.1 <- density(nico_morpholin$x)
plot(sp.1, main = "Kernel Density of S2_x")
polygon(sp.1, col = "red", border = "black")
rug(nico_morpholin$x, col = "brown")

sp.2 <- density(nico_morpholin$y)
plot(sp.2, main = "Kernel Density of S2_y")
polygon(sp.2, col = "red", border = "black")
rug(nico_morpholin$y, col = "brown")

sp.3 <- density(nico_morpholin$z)
plot(sp.3, main = "Kernel Density of S2_z")
polygon(sp.3, col = "red", border = "black")
rug(nico_morpholin$z, col = "brown")

##cola ancho
aw <- density(nico_morpholin$AW)
plot(aw, main = "Kernel Density of AW")
polygon(aw, col = "red", border = "black")
rug(nico_morpholin$AW, col = "brown")

#areola ancho
arew <- density(nico_morpholin$AreW)
plot(arew, main = "Kernel Density of AreW")
polygon(arew, col = "red", border = "black")
rug(nico_morpholin$AreW, col = "brown")

#areola largo
arel <- density(nico_morpholin$AreL)
plot(arel, main = "Kernel Density of AreL")
polygon(arel, col = "red", border = "black")
rug(nico_morpholin$AreL, col = "brown")


###Tests of independence of categorical variables###

library(vcd)

#Clado_vs_K4
mytable_1 <- xtabs(~Clado+K4, data = nico_morpholin)
test1 <- chisq.test(mytable_1)
#	Pearson's Chi-squared test
#data:  mytable_1
#X-squared = 80, df = 6, p-value = 3.573e-15
###Aren't independent, there appears to be a relantionship between Clado and K4
#to find the expected frequencies
test1$expected
test1$observed

#Clado_vs_BFD_gdi
mytable_2 <- xtabs(~Clado+BFD_gdi, data = nico_morpholin)
chisq.test(mytable_2)
###~Region+K7###
##Pearson's Chi-squared test
##data:  mytable_2
##X-squared = 80, df = 12, p-value = 4.127e-12
##Aren't independent, there appears to be a relantionship between Clado and K4
##p < .01. The p-values are the probability of obtaining the sampled results, assuming independence of the row and column variables in the population.

#Clado_vs_COI
mytable_3 <- xtabs(~Clado+COI, data = nico_morpholin)
chisq.test(mytable_3)
###~Region+X4_species###
##	Pearson's Chi-squared test
##data:  mytable_3
##X-squared = 80, df = 10, p-value = 5.02e-13
##Aren't independent, there appears to be a relantionship between clados and COI

#Clado_vs_K5
mytable_4 <- xtabs(~Clado+K5, data = nico_morpholin)
chisq.test(mytable_4)
###~Region+X4_species###
##	Pearson's Chi-squared test
##data:  mytable_4
##X-squared = 80, df = 8, p-value = 4.889e-14
##Aren't independent, there appears to be a relantionship between clados and K5

#K4 vs COI
mytable_5 <- xtabs(~K4+COI, data = nico_morpholin)
chisq.test(mytable_5)
##Pearson's Chi-squared test
##data:  mytable_5
##X-squared =120, df = 18, p-value < 2.2e-16

#K4 vs BFD_gdi
mytable_6 <- xtabs(~K4+BFD_gdi, data = nico_morpholin)
chisq.test(mytable_6)
##Pearson's Chi-squared test
##data:  mytable_7
##X-squared =120, df = 18, p-value < 2.2e-16

#K4 vs K5
mytable_7 <- xtabs(~K4+K5, data = nico_morpholin)
chisq.test(mytable_7)
##Pearson's Chi-squared test
##data:  mytable_7
##X-squared = 120, df = 12, p-value < 2.2e-16

#COI vs BFD_gdi
mytable_8 <- xtabs(~COI+BFD_gdi, data = nico_morpholin) 
chisq.test(mytable_8)
##Pearson's Chi-squared test
##data:  mytable_8
#X-squared = 200, df = 30, p-value < 2.2e-16

#COI vs K5
mytable_9 <- xtabs(~COI+K5, data = nico_morpholin)
chisq.test(mytable_9) #Are not independent
#Pearson's Chi-squared test
#data:  mytable_9
#X-squared = 160, df = 20, p-value < 2.2e-16

#BFD_gdi vs K5
mytable_10 <- xtabs(~BFD_gdi+K5, data = nico_morpholin)
chisq.test(mytable_10) #Are not independent
#Pearson's Chi-squared test
#data:  mytable_10
#X-squared = 160, df = 24, p-value < 2.2e-16


#############################################
###linear models
#########################
modelo1 <- aov(TL+CL+Cel+AW+AreW+AreL+RL+RW+DL+PrL+PrW+x+y+z ~ Clado, data= nico_morpholin)
par(mfrow=c(3,2))
plot(modelo1)#Primero vemos los supuestos. Si todo esta bien, entonces seguimos
hist(resid(modelo1))
boxplot(resid(modelo1)~modelo1$fitted, main="Residuals vs. Fitted")#homogeneidad de varianzas
shapiro.test(modelo1$residuals) 
#Shapiro-Wilk normality test
#data:  modelo1$residuals
#W = 0.98344, p-value = 0.8139
#From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
summary(modelo1)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2   1504   752.0   1.177   0.32
#Residuals   37  23647   639.1  
summary.lm(modelo1)
#Coefficients:
#            Estimate Std. Error   t value Pr(>|t|)    
#(Intercept)   253.53       5.80  43.713   <2e-16 ***
#CladoNorth    -11.63       9.10  -1.278    0.209    
#CladoSouth    -13.20      10.65  -1.239    0.223    
#Residual standard error: 25.28 on 37 degrees of freedom
#Multiple R-squared:  0.0598,	Adjusted R-squared:  0.008976 
#F-statistic: 1.177 on 2 and 37 DF,  p-value: 0.3196
TukeyHSD(modelo1)
#                    diff       lwr      upr     p adj
#North-Central -11.633320 -33.84963 10.58299 0.4160537
#South-Central -13.198224 -39.21213 12.81568 0.4384205
#South-North    -1.564904 -29.30056 26.17075 0.9895958
plot(TukeyHSD(modelo1))

modelo2 <- aov(TL+CL+Cel+AW+AreW+AreL+RL+RW+DL+PrL+PrW+x+y+z ~ K4, data= nico_morpholin)
par(mfrow=c(3,2))
plot(modelo2)#Primero vemos los supuestos. Si todo esta bien, entonces seguimos
hist(resid(modelo2))
boxplot(resid(modelo2)~modelo2$fitted, main="Residuals vs. Fitted")#homogeneidad de varianzas
shapiro.test(modelo2$residuals) 
#Shapiro-Wilk normality test
#data:  modelo2$residuals
#W = 0.99087, p-value = 0.984
#From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
summary(modelo2)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K4           3   3992  1330.6   2.264 0.0977 .
#Residuals   36  21160   587.8 
TukeyHSD(modelo2)
#                     diff       lwr      upr     p adj
#North1-Central -36.889474 -77.45443  3.67548 0.0858009
#North2-Central  -4.056474 -29.56582 21.45287 0.9732498
#South-Central  -13.198224 -40.71751 14.32106 0.5740681
#North2-North1   32.833000 -10.14915 75.81515 0.1866175 ojo
#South-North1    23.691250 -20.51338 67.89588 0.4814158
#South-North2    -9.141750 -40.11368 21.83018 0.8562604

##########modelo 3#######################################
modelo3 <- aov(TL+CL+Cel+AW+AreW+AreL+RL+RW+DL+PrL+PrW+x+y+z ~ K5, data= nico_morpholin)
plot(modelo3)
hist(resid(modelo3))
boxplot(resid(modelo3)~modelo3$fitted, main="Residuals vs. Fitted")
shapiro.test(modelo3$residuals)
#Shapiro-Wilk normality test
#data:  modelo3$residuals
#W = 0.98716, p-value = 0.9232
summary(modelo3)
#             Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   4169  1042.2   1.739  0.164
#Residuals   35  20983   599.5 
summary.lm(modelo3)
#Coefficients:
#           Estimate Std.  Error  t value  Pr(>|t|)    
#(Intercept)  256.749      8.162  31.458   <2e-16 ***
#K5Center2     -6.117     11.250  -0.544   0.5901    
#K5North1     -40.109     16.323  -2.457   0.0191 *  
#K5North2      -7.276     11.250  -0.647   0.5220    
#K5South1     -16.418     11.897  -1.380   0.1764    

#Residual standard error: 24.48 on 35 degrees of freedom
#Multiple R-squared:  0.1658,	Adjusted R-squared:  0.07041 
#F-statistic: 1.739 on 4 and 35 DF,  p-value: 0.1636
 
TukeyHSD(modelo3)
#                      diff       lwr       upr     p adj
#Center2-Center1  -6.116889 -38.46116 26.227380 0.9820042
#North1-Center1  -40.108889 -87.03889  6.821108 0.1241500 ojo
#North2-Center1   -7.275889 -39.62016 25.068380 0.9661169
#South1-Center1  -16.417639 -50.62346 17.788181 0.6441573
#North1-Center2  -33.992000 -80.33166 12.347659 0.2391091
#North2-Center2   -1.159000 -32.64060 30.322599 0.9999701
#South1-Center2  -10.300750 -43.69203 23.090528 0.8998139
#North2-North1    32.833000 -13.50666 79.172659 0.2701348
#South1-North1    23.691250 -23.96639 71.348888 0.6135419
#South1-North2    -9.141750 -42.53303 24.249528 0.9326522

###############modelo 4#########################################
modelo4 <- aov(TL+CL+Cel+AW+AreW+AreL+RL+RW+DL+PrL+PrW+x+y+z ~ COI, data= nico_morpholin)
par(mfrow=c(3,2))
plot(modelo4)
hist(resid(modelo4))
boxplot(resid(modelo4)~modelo4$fitted, main="Residuals vs. Fitted")
shapiro.test(modelo4$residuals)
#Shapiro-Wilk normality test
#data:  modelo4$residuals
#W = 0.96849, p-value = 0.322
summary(modelo4)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5   6717  1343.4   2.478 0.0512 .
#Residuals   34  18435   542.2
summary.lm(modelo4)
#Coefficients:
#           Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  256.749      7.762  33.079   <2e-16 ***
#COICenter2    -6.117     10.699  -0.572   0.5713    
#COINorth1    -40.109     15.523  -2.584   0.0142 *  
#COINorth2     -7.276     10.699  -0.680   0.5011    
#COISouth1    -34.264     13.993  -2.449   0.0196 *  
#COISouth2      1.429     13.993   0.102   0.9193    

#Residual standard error: 23.29 on 34 degrees of freedom
#Multiple R-squared:  0.2671,	Adjusted R-squared:  0.1593 
#F-statistic: 2.478 on 5 and 34 DF,  p-value: 0.05119

TukeyHSD(modelo4)
#                      diff       lwr       upr     p adj
#Center2-Center1  -6.116889 -38.40831 26.174533 0.9922459
#North1-Center1  -40.108889 -86.96221  6.744431 0.1290860 ojo
#North2-Center1   -7.275889 -39.56731 25.015533 0.9830000
#South1-Center1  -34.263889 -76.49690  7.969123 0.1684967 ojo
#South2-Center1    1.428611 -40.80440 43.661623 0.9999983
#North1-Center2  -33.992000 -80.25595 12.271946 0.2563334
#North2-Center2   -1.159000 -32.58916 30.271162 0.9999974
#South1-Center2  -28.147000 -69.72520 13.431197 0.3402719
#South2-Center2    7.545500 -34.03270 49.123697 0.9936380
#North2-North1    32.833000 -13.43095 79.096946 0.2909069
#South1-North1     5.845000 -47.83222 59.522221 0.9994424
#South2-North1    41.537500 -12.13972 95.214721 0.2081414
#South1-North2   -26.988000 -68.56620 14.590197 0.3858115
#South2-North2     8.704500 -32.87370 50.282697 0.9877681
#South2-South1    35.692500 -14.00295 85.387950 0.2788059

###########################modelo 5############################
modelo5 <- aov(TL+CL+Cel+AW+AreW+AreL+RL+RW+DL+PrL+PrW+x+y+z ~ BFD_gdi, data= nico_morpholin)
par(mfrow=c(3,2))
plot(modelo5)
hist(resid(modelo5))
boxplot(resid(modelo5)~modelo5$fitted, main="Residuals vs. Fitted")
shapiro.test(modelo5$residuals)
#Shapiro-Wilk normality test
#data:  modelo5$residuals
#W = 0.96988, p-value = 0.3567
summary(modelo5)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#BFD_gdi      6   6866  1144.3   2.065 0.0845 .
#Residuals   33  18286   554.1 
summary.lm(modelo5)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     256.749      7.847  32.721   <2e-16 ***
#BFD_gdiCenter2   -6.117     10.816  -0.566   0.5755    
#BFD_gdiNorth1   -40.109     15.693  -2.556   0.0154 *  
#BFD_gdiNorth2    -7.276     10.816  -0.673   0.5058    
#BFD_gdiSouth1   -34.264     14.145  -2.422   0.0211 *  
#BFD_gdiSouth2    12.001     24.813   0.484   0.6318    
#BFD_gdiSouth3    -2.096     15.693  -0.134   0.8946    

#Residual standard error: 23.54 on 33 degrees of freedom
#Multiple R-squared:  0.273,	Adjusted R-squared:  0.1408 
#F-statistic: 2.065 on 6 and 33 DF,  p-value: 0.08447

TukeyHSD(modelo5)
#                      diff       lwr       upr     p adj
#Center2-Center1  -6.116889 -40.04632  27.812547 0.9973836
#North1-Center1  -40.108889 -89.33889   9.121111 0.1728665 ojo
#North2-Center1   -7.275889 -41.20532  26.653547 0.9932696
#South1-Center1  -34.263889 -78.63921  10.111433 0.2215558
#South2-Center1   12.001111 -65.83835  89.840576 0.9989098
#South3-Center1   -2.095556 -51.32556  47.134444 0.9999994
#North1-Center2  -33.992000 -82.60273  14.618730 0.3259829
#North2-Center2   -1.159000 -34.18349  31.865488 0.9999998
#South1-Center2  -28.147000 -71.83429  15.540291 0.4211800
#South2-Center2   18.118000 -59.33129  95.567289 0.9893338
#South3-Center2    4.021333 -44.58940  52.632063 0.9999705
#North2-North1    32.833000 -15.77773  81.443730 0.3657086
#South1-North1     5.845000 -50.55505  62.245050 0.9998894
#South2-North1    52.110000 -33.15886 137.378861 0.4838646
#South3-North1    38.013333 -22.28086  98.307523 0.4469143
#South1-North2   -26.988000 -70.67529  16.699291 0.4710902
#South2-North2    19.277000 -58.17229  96.726289 0.9852875
#South3-North2     5.180333 -43.43040  53.791063 0.9998698
#South2-South1    46.265000 -36.29622 128.826220 0.5840596
#South3-South1    32.168333 -24.23172  88.568384 0.5641555
#South3-South2   -14.096667 -99.36553  71.172194 0.9983852


##################
####Largo Total###
##################
fit1 <- aov(TL+CL+Cel ~ BFD_gdi, data= nico_morpholin)
par(mfrow=c(3,2))
plot(fit1)
hist(resid(fit1))
boxplot(resid(fit1)~fit1$fitted, main="Residuals vs. Fitted")
shapiro.test(fit1$residuals)
#	Shapiro-Wilk normality test
#data:  fit1$residuals
#W = 0.93864, p-value = 0.03111
fit1a <- aov(TL ~ BFD_gdi, data= nico_morpholin)
fit1b <- aov(CL ~ BFD_gdi, data= nico_morpholin)
fit1c <- aov(Cel ~ BFD_gdi, data= nico_morpholin)
summary(fit1c)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#BFD_gdi      6  67.04  11.174    2.12 0.0773 .
#Residuals   33 173.94   5.271 
TukeyHSD(fit1c)
summary(fit1b)
#             Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  137.9   22.98    1.92  0.107
#Residuals   33  395.0   11.97
summary(fit1a)
#          Df Sum Sq Mean Sq   F value Pr(>F)  
#BFD_gdi      6    581   96.84   2.214 0.0664 .
#Residuals   33   1444   43.74  
summary(fit1)
##            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6   1922   320.4   2.137  0.0752
#Residuals   33   4948   149.9               
#
#Residual standard error: 12.24 on 33 degrees of freedom
#Multiple R-squared:  0.2798,	Adjusted R-squared:  0.1489 
#F-statistic: 2.137 on 6 and 33 DF,  p-value: 0.07522
#4 observations deleted due to missingness
###The ANOVA F test for Region is not significant (p > .05), providing evidence that 
###in candidate species (K = 7) doesn't exists differences in total lenght
summary.lm(fit1)
TukeyHSD(fit1)
plot(TukeyHSD(fit1))
#                       diff       lwr       upr     p adj
#Center2-Center1  -4.8360000 -22.48499 12.812992 0.9760767
#North1-Center1  -19.8200000 -45.42785  5.787849 0.2192611
#North2-Center1   -4.8740000 -22.52299 12.774992 0.9751288
#South1-Center1  -20.4575000 -43.54010  2.625103 0.1101896 ojo
#South2-Center1    2.5700000 -37.91956 43.059565 0.9999939
#South3-Center1   -0.5333333 -26.14118 25.074516 1.0000000
#North1-Center2  -14.9840000 -40.26973 10.301725 0.5200850
#North2-Center2   -0.0380000 -17.21627 17.140267 1.0000000
#South1-Center2  -15.6215000 -38.34621  7.103212 0.3454019
#South2-Center2    7.4060000 -32.88061 47.692608 0.9970861
#South3-Center2    4.3026667 -20.98306 29.588392 0.9981022
#North2-North1    14.9460000 -10.33973 40.231725 0.5230472
#South1-North1    -0.6375000 -29.97498 28.699977 1.0000000
#South2-North1    22.3900000 -21.96410 66.744096 0.6932041
#South3-North1    19.2866667 -12.07642 50.649749 0.4764851
#South1-North2   -15.5835000 -38.30821  7.141212 0.3481945
#South2-North2     7.4440000 -32.84261 47.730608 0.9970027
#South3-North2     4.3406667 -20.94506 29.626392 0.9980068
#South2-South1    23.0275000 -19.91817 65.973169 0.6320957
#South3-South1    19.9241667  -9.41331 49.261644 0.3593615
#South3-South2    -3.1033333 -47.45743 41.250763 0.9999891

fit1.1 <- aov(TL+CL+Cel ~ COI, data= nico_morpholin)
fit1.1A <- aov(TL ~ COI, data= nico_morpholin)
fit1.1b <- aov(CL ~ COI, data= nico_morpholin)
fit1.1c <- aov(Cel ~ COI, data= nico_morpholin)
summary(fit1.1c)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5  66.51  13.303   2.592 0.0432 *
#Residuals   34 174.47   5.131 
summary(fit1.1b)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5  134.2   26.85    2.29 0.0675 .
#Residuals   34  398.6   11.72 
summary(fit1.1A)
#           Df  Sum Sq Mean Sq F value  Pr(>F)  
#COI          5    581  116.20   2.737  0.035 *
#Residuals   34   1444   42.46 
summary(fit1.1)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5   1915   383.0   2.628  0.041 *
#Residuals   34   4955   145.7  
summary.lm(fit1.1)
#Residual standard error: 12.07 on 34 degrees of freedom
#Multiple R-squared:  0.2788,	Adjusted R-squared:  0.1727 
#F-statistic: 2.628 on 5 and 34 DF,  p-value: 0.04103
TukeyHSD(fit1.1)
plot(TukeyHSD(fit1.1))
#                     diff        lwr       upr     p adj
#Center2-Center1  -4.8360 -21.577121 11.905121 0.9506641
#North1-Center1  -19.8200 -44.110571  4.470571 0.1640110
#North2-Center1   -4.8740 -21.615121 11.867121 0.9490530
#South1-Center1  -20.4575 -42.352725  1.437725 0.0783117
#South2-Center1    0.2425 -21.652725 22.137725 1.0000000
#North1-Center2  -14.9840 -38.969017  9.001017 0.4279327
#North2-Center2   -0.0380 -16.332611 16.256611 1.0000000
#South1-Center2  -15.6215 -37.177244  5.934244 0.2698433
#South2-Center2    5.0785 -16.477244 26.634244 0.9792976
#North2-North1    14.9460  -9.039017 38.931017 0.4307377
#South1-North1    -0.6375 -28.465845 27.190845 0.9999998
#South2-North1    20.0625  -7.765845 47.890845 0.2750291
#South1-North2   -15.5835 -37.139244  5.972244 0.2722635
#South2-North2     5.1165 -16.439244 26.672244 0.9786093
#South2-South1    20.7000  -5.064042 46.464042 0.1762958

fit1.2 <- aov(TL+CL+Cel ~ K4, data= nico_morpholin)
fit1.2a <- aov(TL ~ K4, data= nico_morpholin)
fit1.2b <- aov(CL ~ K4, data= nico_morpholin)
fit1.2c <- aov(Cel ~ K4, data= nico_morpholin)
summary(fit1.2c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  31.22  10.406   1.786  0.167
#Residuals   36 209.76   5.827
summary(fit1.2b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3   63.9   21.30   1.635  0.198
#Residuals   36  469.0   13.03
summary(fit1.2a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  298.9   99.64   2.079   0.12
#Residuals   36 1725.6   47.93
summary(fit1.2)
#             Df Sum Sq Mean Sq F value Pr(>F)
#K4           3    947   315.7   1.919  0.144
#Residuals   36   5923   164.5 
plot(fit1.2)
summary.lm(fit1.2)
#Residual standard error: 12.83 on 36 degrees of freedom
#Multiple R-squared:  0.1379,	Adjusted R-squared:  0.06604 
#F-statistic: 1.919 on 3 and 36 DF,  p-value: 0.1439
TukeyHSD(fit1.2)
plot(TukeyHSD(fit1.2))
#                    diff        lwr       upr     p adj
#North1-Central -17.274737 -38.735812  4.186338 0.1517224
#North2-Central  -2.328737 -15.824572 11.167098 0.9662615
#South-Central   -7.562237 -22.121442  6.996968 0.5082237
#North2-North1   14.946000  -7.793903 37.685903 0.3039995
#South-North1     9.712500 -13.674165 33.099165 0.6806681
#South-North2    -5.233500 -21.619344 11.152344 0.8251254

fit1.3 <- aov(TL+CL+Cel ~ Clado, data= nico_morpholin)
fit1.3a <- aov(TL ~ Clado, data= nico_morpholin)
fit1.3b <- aov(CL ~ Clado, data= nico_morpholin)
fit1.3c <- aov(Cel ~ Clado, data= nico_morpholin)
summary(fit1.3c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2  14.91   7.454    1.22  0.307
#Residuals   37 226.08   6.110
summary(fit1.3b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2   24.1   12.03   0.875  0.425
#Residuals   37  508.8   13.75
summary(fit1.3a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2  146.3   73.16   1.441   0.25
#Residuals   37 1878.2   50.76  
summary(fit1.3)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#Region       2    432   215.9   1.241  0.301
#Residuals   37   6438   174.0  
TukeyHSD(fit1.3)
plot(TukeyHSD(fit1.3))
#                   diff       lwr       upr     p adj
#North-Central -5.777814 -17.36981  5.814182 0.4508932
#South-Central -7.562237 -21.13574  6.011262 0.3718489
#South-North   -1.784423 -16.25629 12.687448 0.9513364

fit1.4 <- aov(TL+CL+Cel ~ K5, data= nico_morpholin)
plot(fit1.4)
hist(resid(fit1.4))
boxplot(resid(fit1.4)~fit1.4$fitted, main="Residuals vs. Fitted")
shapiro.test(fit1.4$residuals)
#Shapiro-Wilk normality test
#data:  fit1.4$residuals
#W = 0.97075, p-value = 0.3802
fit1.4a <- aov(TL ~ K5, data= nico_morpholin)
fit1.4b <- aov(CL ~ K5, data= nico_morpholin)
fit1.4c <- aov(Cel ~ K5, data= nico_morpholin)
summary(fit1.4c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4  34.31   8.578   1.453  0.237
#Residuals   35 206.67   5.905  
summary(fit1.4b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   70.9   17.73   1.343  0.274
#Residuals   35  462.0   13.20
summary(fit1.4a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4  336.4   84.09   1.743  0.163
#Residuals   35 1688.2   48.23
summary(fit1.4)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   1058   264.5   1.593  0.198
#Residuals   35   5812   166.1 
summary.lm(fit1.4)
#Residual standard error: 12.89 on 35 degrees of freedom
#Multiple R-squared:  0.154,	Adjusted R-squared:  0.05732 
#F-statistic: 1.593 on 4 and 35 DF,  p-value: 0.1979
TukeyHSD(fit1.4)
#                    diff        lwr       upr     p adj
#Center2-Center1  -4.8360 -21.858533 12.186533 0.9236887
#North1-Center1  -19.8200 -44.518887  4.878887 0.1667497
#North2-Center1   -4.8740 -21.896533 12.148533 0.9216605
#South1-Center1  -10.1075 -28.109753  7.894753 0.4986077
#North1-Center2  -14.9840 -39.372197  9.404197 0.4085817
#North2-Center2   -0.0380 -16.606517 16.530517 1.0000000
#South1-Center2   -5.2715 -22.845066 12.302066 0.9085794
#North2-North1    14.9460  -9.442197 39.334197 0.4111298
#South1-North1     9.7125 -15.369338 34.794338 0.7984417
#South1-North2    -5.2335 -22.807066 12.340066 0.9107316


########################################################
###analizando la pleura 2
S2 <- c("x", "y", "z")
head (nico_morpholin[S2])
summary(nico_morpholin[S2])

fit2 <- aov(x+y+z ~ BFD_gdi, data= nico_morpholin)
plot(fit2)
hist(resid(fit2))
boxplot(resid(fit2)~fit2$fitted, main="Residuals vs. Fitted")
shapiro.test(fit2$residuals)
#	Shapiro-Wilk normality test
#data:  fit2$residuals
#W = 0.98299, p-value = 0.7981
fit2a <- aov(x ~ BFD_gdi, data= nico_morpholin)
fit2b <- aov(y ~ BFD_gdi, data= nico_morpholin)
fit2c <- aov(z ~ BFD_gdi, data= nico_morpholin)
summary(fit2c)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#BFD_gdi      6  8.654  1.4423   2.312 0.0566 .
#Residuals   33 20.584  0.6238 
TukeyHSD(fit2c)
summary(fit2b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  9.702  1.6169   1.825  0.125
#Residuals   33 29.245  0.8862 
summary(fit2a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  6.473  1.0788   1.927  0.106
#Residuals   33 18.471  0.5597
summary(fit2)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi          6  63.88  10.646   2.089  0.0813 .
#Residuals   33 168.18   5.096  
summary.lm(fit2)
#Residual standard error: 2.257 on 33 degrees of freedom
#Multiple R-squared:  0.2753,	Adjusted R-squared:  0.1435 
#F-statistic: 2.089 on 6 and 33 DF,  p-value: 0.08127
TukeyHSD(fit2)
#                       diff         lwr        upr     p adj
#Center2-Center1 -0.56411111  -3.8180167  2.6897945 0.9978944
#North1-Center1  -4.02777778  -8.7490402  0.6934847 0.1364344 ojo
#North2-Center1  -0.02511111  -3.2790167  3.2287945 1.0000000
#South1-Center1  -2.66361111  -6.9192996  1.5920774 0.4555549
#South2-Center1   1.76888889  -5.6960825  9.2338603 0.9885935
#South3-Center1  -0.66111111  -5.3823736  4.0601514 0.9993693
#North1-Center2  -3.46366667  -8.1255398  1.1982065 0.2601542
#North2-Center2   0.53900000  -2.6281191  3.7061191 0.9981007
#South1-Center2  -2.09950000  -6.2892048  2.0902048 0.7002240
#South2-Center2   2.33300000  -5.0945528  9.7605528 0.9537100
#South3-Center2  -0.09700000  -4.7588731  4.5648731 1.0000000
#North2-North1    4.00266667  -0.6592065  8.6645398 0.1317754 ojo
#South1-North1    1.36416667  -4.0447190  6.7730523 0.9842562
#South2-North1    5.79666667  -2.3807998 13.9741331 0.3107431
#South3-North1    3.36666667  -2.4156753  9.1490087 0.5405557
#South1-North2   -2.63850000  -6.8282048  1.5512048 0.4482458
#South2-North2    1.79400000  -5.6335528  9.2215528 0.9874010
#South3-North2   -0.63600000  -5.2978731  4.0258731 0.9994566
#South2-South1    4.43250000  -3.4852979 12.3502979 0.5851704
#South3-South1    2.00250000  -3.4063857  7.4113857 0.9033724
#South3-South2   -2.43000000 -10.6074665  5.7474665 0.9644252


fit2.1 <- aov(x+y+z ~ COI, data= nico_morpholin)
fit2.1a <- aov(x ~ COI, data= nico_morpholin)
fit2.1b <- aov(y ~ COI, data= nico_morpholin)
fit2.1c <- aov(z ~ COI, data= nico_morpholin)
summary(fit2.1c)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5  8.327   1.665   2.708 0.0365 *
#Residuals   34 20.911   0.615 
TukeyHSD(fit2.1c)
summary(fit2.1b)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5  9.698  1.9396   2.255 0.0711 .
#Residuals   34 29.249  0.8603
TukeyHSD(fit2.1b)
summary(fit2.1a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5  4.305   0.861   1.418  0.243
#Residuals   34 20.639   0.607 
summary(fit2.1)
#             Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5  59.45  11.890   2.342 0.0625 .
#Residuals   34 172.61   5.077  
summary.lm(fit2.1)
#Residual standard error: 2.253 on 34 degrees of freedom
#Multiple R-squared:  0.2562,	Adjusted R-squared:  0.1468 
#F-statistic: 2.342 on 5 and 34 DF,  p-value: 0.06252
TukeyHSD(fit2.1)
#                     diff         lwr        upr     p adj
#Center2-Center1 -0.56411111 -3.6887299 2.5605077 0.9937891
#North1-Center1  -4.02777778 -8.5614513 0.5058957 0.1055129
#North2-Center1  -0.02511111 -3.1497299 3.0995077 1.0000000
#South1-Center1  -2.66361111 -6.7502092 1.4229869 0.3812679
#South2-Center1  -0.05361111 -4.1402092 4.0329869 1.0000000
#North1-Center2  -3.46366667 -7.9403105 1.0129772 0.2082774
#North2-Center2   0.53900000 -2.5022806 3.5802806 0.9943002
#South1-Center2  -2.09950000 -6.1227361 1.9237361 0.6198396
#South2-Center2   0.51050000 -3.5127361 4.5337361 0.9988307
#North2-North1    4.00266667 -0.4739772 8.4793105 0.1017503
#South1-North1    1.36416667 -3.8298088 6.5581422 0.9668742
#South2-North1    3.97416667 -1.2198088 9.1681422 0.2182429
#South1-North2   -2.63850000 -6.6617361 1.3847361 0.3745509
#South2-North2   -0.02850000 -4.0517361 3.9947361 1.0000000
#South2-South1    2.61000000 -2.1986869 7.4186869 0.5800665


fit2.2 <-  aov(x+y+z ~ K4, data= nico_morpholin)
fit2.2a <- aov(x ~ K4, data= nico_morpholin)
fit2.2b <- aov(y ~ K4, data= nico_morpholin)
fit2.2c <- aov(z ~ K4, data= nico_morpholin)
summary(fit2.2c)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K4           3  6.697  2.2322   3.565 0.0234 *
#Residuals   36 22.541  0.6262
TukeyHSD(fit2.2c)
#                       diff         lwr         upr     p adj
#North1-Central -1.419649123 -2.74364500 -0.09565324 0.0316468*
#North2-Central -0.004315789 -0.83691294  0.82828136 0.9999990
#South-Central  -0.572565789 -1.47076529  0.32563371 0.3300604
#North2-North1   1.415333333  0.01244285  2.81822382 0.0473275*
#South-North1    0.847083333 -0.59570773  2.28987439 0.4017842
#South-North2   -0.568250000 -1.57914018  0.44264018 0.4399471
summary(fit2.2b)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K4           3  7.381  2.4604   2.806 0.0534 .
#Residuals   36 31.566  0.8768 
TukeyHSD(fit2.2b)
#                    diff         lwr         upr     p adj
#North1-Central -1.60649123 -3.17325783 -0.03972463 0.0426942*
#North2-Central  0.01884211 -0.96642184  1.00410605 0.9999500
#South-Central  -0.36315789 -1.42605318  0.69973739 0.7942504
#North2-North1   1.62533333 -0.03479415  3.28546082 0.0568496 .
#South-North1    1.24333333 -0.46401098  2.95067765 0.2214979
#South-North2   -0.38200000 -1.57824917  0.81424917 0.8252059
summary(fit2.2a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  2.271  0.7571   1.202  0.323
#Residuals   36 22.672  0.6298 
summary(fit2.2)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K4           3  44.32  14.772   2.833 0.0519 .
#Residuals   36 187.74   5.215 
summary.lm(fit2.2)
#Residual standard error: 2.284 on 36 degrees of freedom
#Multiple R-squared:  0.191,	Adjusted R-squared:  0.1236 
#F-statistic: 2.833 on 3 and 36 DF,  p-value: 0.05188
TukeyHSD(fit2.2)
#                     diff       lwr       upr     p adj
#North1-Central -3.7308772 -7.55181757 0.09006319 0.0577584
#North2-Central  0.2717895 -2.13101564 2.67459458 0.9900170
#South-Central  -1.0617105 -3.65383849 1.53041743 0.6899389
#North2-North1   4.0026667 -0.04595689 8.05129022 0.0536138
#South-North1    2.6691667 -1.49460658 6.83293991 0.3252462
#South-North2   -1.3335000 -4.25084377 1.58384377 0.6116185
###ojo con las relaciones entre North1 y Central, North1 vs North2


fit2.3 <- aov(x+y+z ~ K5, data= nico_morpholin)
fit2.3a <- aov(x ~ K5, data= nico_morpholin)
fit2.3b <- aov(y ~ K5, data= nico_morpholin)
fit2.3c <- aov(z ~ K5, data= nico_morpholin)
summary(fit2.3c)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K5           4  6.698   1.675     2.6 0.0528 .
#Residuals   35 22.540   0.644
TukeyHSD(fit2.3c)
summary(fit2.3b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4  7.387  1.8467   2.048  0.109
#Residuals   35 31.560  0.9017
summary(fit2.3a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4  3.505  0.8762    1.43  0.244
#Residuals   35 21.439  0.6125
summary(fit2.3)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K5           4  45.82  11.456   2.153 0.0948 .
#Residuals   35 186.23   5.321 
summary.lm(fit2.3)
#Residual standard error: 2.307 on 35 degrees of freedom
#Multiple R-squared:  0.1975,	Adjusted R-squared:  0.1058 
#F-statistic: 2.153 on 4 and 35 DF,  p-value: 0.09484
TukeyHSD(fit2.3)
#                       diff        lwr       upr     p adj
#Center2-Center1 -0.56411111 -3.6112522 2.4830300 0.9833726
#North1-Center1  -4.02777778 -8.4490347 0.3934792 0.0885333 .
#North2-Center1  -0.02511111 -3.0722522 3.0220300 0.9999999
#South1-Center1  -1.35861111 -4.5811282 1.8639060 0.7444465
#North1-Center2  -3.46366667 -7.8293081 0.9019748 0.1751937
#North2-Center2   0.53900000 -2.4268693 3.5048693 0.9844784
#South1-Center2  -0.79450000 -3.9402795 2.3512795 0.9489789
#North2-North1    4.00266667 -0.3629748 8.3683081 0.0853839 .
#South1-North1    2.66916667 -1.8206410 7.1589743 0.4416380
#South1-North2   -1.33350000 -4.4792795 1.8122795 0.7406667


fit2.4 <- aov(x+y+z ~ Clado, data= nico_morpholin)
fit2.4a <- aov(x ~ Clado, data= nico_morpholin)
fit2.4b <- aov(y ~ Clado, data= nico_morpholin)
fit2.4c <- aov(z ~ Clado, data= nico_morpholin)
summary(fit2.4c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2  2.074  1.0369   1.412  0.256
#Residuals   37 27.164  0.7342
summary(fit2.4b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2   1.28  0.6425   0.631  0.538
#Residuals   37  37.66  1.0179 
summary(fit2.4a)
#Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2  0.136  0.0678   0.101  0.904
#Residuals   37 24.808  0.6705
summary(fit2.4)
#           Df Sum Sq Mean Sq F value Pr(>F)
#Clado       2   7.34   3.672   0.605  0.552
#Residuals   37 224.71   6.073   
TukeyHSD(fit2.4)
#no hay diferencias
#                  diff       lwr      upr     p adj
#North-Central -0.6519028 -2.817562 1.513756 0.7444791
#South-Central -1.0617105 -3.597561 1.474140 0.5677634
#South-North   -0.4098077 -3.113495 2.293880 0.9274403


#########################################################
###Evaluando la quela#########
##########################################
fit3 <- aov(DL+PrL+PrW ~ BFD_gdi, data= nico_morpholin)
fit3a <- aov(DL ~ BFD_gdi, data= nico_morpholin)
fit3b <- aov(PrL ~ BFD_gdi, data= nico_morpholin)
fit3c <- aov(PrW ~ BFD_gdi, data= nico_morpholin)
summary(fit3c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  31.81   5.301   1.383  0.251
#Residuals   33 126.50   3.833
summary(fit3b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  113.6   18.93   1.522  0.202
#Residuals   33  410.5   12.44 
summary(fit3a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  52.86   8.810   1.875  0.115
#Residuals   33 155.07   4.699
summary(fit3)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K_7          6  522.1   87.02   1.612  0.175
#Residuals   33 1781.9   54.00  
summary.lm(fit3)
#Residual standard error: 7.348 on 33 degrees of freedom
#Multiple R-squared:  0.2266,	Adjusted R-squared:  0.08599 
#F-statistic: 1.612 on 6 and 33 DF,  p-value: 0.1751
TukeyHSD(fit3)
#                      diff        lwr       upr     p adj
#Center2-Center1  -1.2837778 -11.875497  9.307942 0.9997250
#North1-Center1  -11.4811111 -26.849193  3.886971 0.2544700
#North2-Center1   -1.1127778 -11.704497  9.478942 0.9998802
#South1-Center1   -6.9777778 -20.830379  6.874824 0.6952718
#South2-Center1    7.8422222 -16.456849 32.141293 0.9474825
#South3-Center1   -0.4244444 -15.792526 14.943637 1.0000000
#North1-Center2  -10.1973333 -25.372098  4.977432 0.3716359
#North2-Center2    0.1710000 -10.138223 10.480223 1.0000000
#South1-Center2   -5.6940000 -19.331820  7.943820 0.8425260
#South2-Center2    9.1260000 -15.051270 33.303270 0.8951424
#South3-Center2    0.8593333 -14.315432 16.034098 0.9999969
#North2-North1    10.3683333  -4.806432 25.543098 0.3523754
#South1-North1     4.5033333 -13.103016 22.109683 0.9830846
#South2-North1    19.3233333  -7.294965 45.941632 0.2846921
#South3-North1    11.0566667  -7.765313 29.878646 0.5302563
#South1-North2    -5.8650000 -19.502820  7.772820 0.8237168
#South2-North2     8.9550000 -15.222270 33.132270 0.9031892
#South3-North2     0.6883333 -14.486432 15.863098 0.9999992
#South2-South1    14.8200000 -10.953057 40.593057 0.5548824
#South3-South1     6.5533333 -11.053016 24.159683 0.9011494
#South3-South2    -8.2666667 -34.884965 18.351632 0.9561090


fit3.1 <- aov(DL+PrL+PrW ~ COI, data= nico_morpholin)
fit3.1a <- aov(DL ~ COI, data= nico_morpholin)
fit3.1b <- aov(PrL ~ COI, data= nico_morpholin)
fit3.1c <- aov(PrW ~ COI, data= nico_morpholin)
summary(fit3.1c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5  26.23   5.247   1.351  0.267
#Residuals   34 132.08   3.885
summary(fit3.1b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5  102.9   20.57    1.66  0.171
#Residuals   34  421.2   12.39
summary(fit3.1a)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5  50.54  10.109   2.184  0.079 .
#Residuals   34 157.38   4.629 
TukeyHSD(fit3.1a)
summary(fit3.1)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5  470.9   94.17   1.747  0.151
#Residuals   34 1833.2   53.92  
summary.lm(fit3.1)
#Residual standard error: 7.343 on 34 degrees of freedom
#Multiple R-squared:  0.2044,	Adjusted R-squared:  0.08735 
#F-statistic: 1.747 on 5 and 34 DF,  p-value: 0.1506
TukeyHSD(fit3.1)
plot(TukeyHSD(fit4.1))
#                     diff        lwr       upr     p adj
#Center2-Center1  -1.283778 -11.466689  8.899134 0.9988664
#North1-Center1  -11.481111 -26.256032  3.293810 0.2044757
#North2-Center1   -1.112778 -11.295689  9.070134 0.9994326
#South1-Center1   -6.977778 -20.295712  6.340156 0.6158683
#South2-Center1    1.642222 -11.675712 14.960156 0.9989813
#North1-Center2  -10.197333 -24.786399  4.391732 0.3065625
#North2-Center2    0.171000  -9.740318 10.082318 0.9999999
#South1-Center2   -5.694000 -18.805442  7.417442 0.7771370
#South2-Center2    2.926000 -10.185442 16.037442 0.9837113
#North2-North1    10.368333  -4.220732 24.957399 0.2894675
#South1-North1     4.503333 -12.423465 21.430132 0.9650129
#South2-North1    13.123333  -3.803465 30.050132 0.2064773
#South1-North2    -5.865000 -18.976442  7.246442 0.7553538
#South2-North2     2.755000 -10.356442 15.866442 0.9875639
#South2-South1     8.620000  -7.051170 24.291170 0.5661993


fit3.2 <- aov(DL+PrL+PrW ~ K4, data= nico_morpholin)
fit3.2a <- aov(DL ~ K4, data= nico_morpholin)
fit3.2b <- aov(PrL ~ K4, data= nico_morpholin)
fit3.2c <- aov(PrW ~ K4, data= nico_morpholin)
summary(fit3.2c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  22.66   7.553   2.004  0.131
#Residuals   36 135.65   3.768
summary(fit3.2b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3   70.5    23.5   1.865  0.153
#Residuals   36  453.6    12.6 
summary(fit3.2a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3   22.7   7.565    1.47  0.239
#Residuals   36  185.2   5.145
summary(fit3.2)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  314.4  104.81   1.896  0.148
#Residuals   36 1989.6   55.27
TukeyHSD(fit3.2) 
#                     diff        lwr       upr     p adj
#North1-Central -10.8054386 -23.244206  1.633329 0.1079017
#North2-Central  -0.4371053  -8.259246  7.385036 0.9987618
#South-Central   -1.9921053 -10.430572  6.446361 0.9197524
#North2-North1   10.3683333  -2.811639 23.548306 0.1665615
#South-North1     8.8133333  -4.741500 22.368166 0.3131453
#South-North2    -1.5550000 -11.052181  7.942181 0.9709277


fit3.3 <- aov(DL+PrL+PrW ~ Clado, data= nico_morpholin)
fit3.3a <- aov(DL ~ Clado, data= nico_morpholin)
fit3.3b <- aov(PrL ~ Clado, data= nico_morpholin)
fit3.3c <- aov(PrW ~ Clado, data= nico_morpholin)
summary(fit3.3c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2    4.7   2.352   0.566  0.572
#Residuals   37  153.6   4.152
summary(fit3.3b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2   10.6   5.311   0.383  0.685
#Residuals   37  513.5  13.878  
summary(fit3.3a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2   8.44   4.221   0.783  0.465
#Residuals   37 199.48   5.391
summary(fit3.3)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Region       2   66.4   33.18   0.549  0.582
#Residuals   37 2237.7   60.48  
TukeyHSD(fit3.3) 
#                  diff       lwr      upr     p adj
#North-Central -2.8297976 -9.663844 4.004249 0.5747205
#South-Central -1.9921053 -9.994344 6.010133 0.8167756
#South-North    0.8376923 -7.694180 9.369565 0.9688457


fit3.4 <- aov(DL+PrL+PrW ~ K5, data= nico_morpholin)
fit3.4a <- aov(DL ~ K5, data= nico_morpholin)
fit3.4b <- aov(PrL ~ K5, data= nico_morpholin)
fit3.4c <- aov(PrW ~ K5, data= nico_morpholin)
summary(fit3.4c)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4  23.07   5.767   1.493  0.226
#Residuals   35 135.24   3.864
summary(fit3.4b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   73.9   18.48   1.437  0.242
#Residuals   35  450.2   12.86  
summary(fit3.4a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   25.2   6.299   1.207  0.325
#Residuals   35  182.7   5.221 
summary(fit3.4)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4  322.2   80.56   1.423  0.247
#Residuals   35 1981.8   56.62 


###################################
###evaluando el rostro########
fit4 <- aov(RL+RW ~ BFD_gdi, data= nico_morpholin)
fit4a <- aov(RL ~ BFD_gdi, data= nico_morpholin)
fit4b <- aov(RW ~ BFD_gdi, data= nico_morpholin)
summary(fit4b)
#            Df Sum Sq Mean Sq F value  Pr(>F)   
#BFD_gdi      6  4.590  0.7650   3.586 0.00759 **
#Residuals   33  7.041  0.2134 
TukeyHSD(fit4b)
#                       diff        lwr         upr     p adj
#North1-Center2  -1.06866667 -2.0225199 -0.11481346 0.0199921 *
#North2-Center2  -0.73600000 -1.3840156 -0.08798435 0.0176913 *
summary(fit4a)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#BFD_gdi      6  6.356  1.0594   3.383 0.0104 *
#Residuals   33 10.333  0.3131 
TukeyHSD(fit4a)
#                      diff         lwr        upr     p adj
#Center2-Center1  0.8850000  0.07842489  1.6915751 0.0239745 *
#North2-Center2  -0.9700000 -1.75506256 -0.1849374 0.0078365 **
summary(fit4)
#             Df Sum Sq Mean Sq F value  Pr(>F)   
#K_7          6  18.81  3.1354   3.478 0.00896 **
#Residuals   33  29.75  0.9016 
TukeyHSD(fit4) 
#                       diff        lwr         upr     p adj
#Center2-Center1  1.11922222 -0.2494219  2.48786635 0.1697357
#North1-Center1  -0.81444444 -2.8002821  1.17139317 0.8531076
#North2-Center1  -0.58677778 -1.9554219  0.78186635 0.8257535
#South1-Center1  -0.37777778 -2.1677876  1.41223206 0.9938210
#South2-Center1  -0.04777778 -3.1876627  3.09210719 1.0000000
#South3-Center1   0.37888889 -1.6069487  2.36472651 0.9964254
#North1-Center2  -1.93366667 -3.8945242  0.02719087 0.0552965 .
#North2-Center2  -1.70600000 -3.0381404 -0.37385963 0.0053651 **
#South1-Center2  -1.49700000 -3.2592561  0.26525607 0.1395968
#South2-Center2  -1.16700000 -4.2911461  1.95714610 0.8996493
#South3-Center2  -0.74033333 -2.7011909  1.22052420 0.8950330
#North2-North1    0.22766667 -1.7331909  2.18852420 0.9997856
#South2-North1    0.76666667 -2.6729050  4.20623832 0.9917309
#South1-North1    0.43666667 -1.8383961  2.71172947 0.9963068
#South3-North1    1.19333333 -1.2388111  3.62547777 0.7198614
#South1-North2    0.20900000 -1.5532561  1.97125607 0.9997575
#South2-North2    0.53900000 -2.5851461  3.66314610 0.9979504
#South3-North2    0.96566667 -0.9951909  2.92652420 0.7164705
#South2-South1    0.33000000 -3.0003509  3.66035093 0.9999149
#South3-South1    0.75666667 -1.5183961  3.03172947 0.9397201
#South3-South2    0.42666667 -3.0129050  3.86623832 0.9996856
plot(TukeyHSD(fit4))


fit4.1 <- aov(RL+RW ~ COI, data= nico_morpholin)
summary(fit4.1)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#COI          5  18.68   3.735   4.249 0.00415 **
#Residuals   34  29.89   0.879 
TukeyHSD(fit4.1)
#                     diff       lwr        upr     p adj
#Center2-Center1  1.1192222 -0.1810448  2.41948927 0.1254087
#North1-Center1  -0.8144444 -2.7010702  1.07218127 0.7813573
#North2-Center1  -0.5867778 -1.8870448  0.71348927 0.7485984
#South1-Center1  -0.3777778 -2.0783592  1.32280366 0.9840395
#South2-Center1   0.2722222 -1.4283592  1.97280366 0.9964602
#North1-Center2  -1.9336667 -3.7965603 -0.07077304 0.0381495 *
#North2-Center2  -1.7060000 -2.9715870 -0.44041299 0.0033448 **
#South1-Center2  -1.4970000 -3.1712142  0.17721424 0.1017317
#South2-Center2  -0.8470000 -2.5212142  0.82721424 0.6498873
#South1-North1    0.4366667 -1.7247346  2.59806796 0.9895902
#North2-North1    0.2276667 -1.6352270  2.09056030 0.9990244
#South2-North1    1.0866667 -1.0747346  3.24806796 0.6557816
#South1-North2    0.2090000 -1.4652142  1.88321424 0.9989191
#South2-North2    0.8590000 -0.8152142  2.53321424 0.6364117
#South2-South1    0.6500000 -1.3510688  2.65106876 0.9210250
plot(TukeyHSD(fit4.1))


fit4.2 <- aov(RL+RW ~ K4, data= nico_morpholin)
fit4.2a <- aov(RL ~ K4, data= nico_morpholin)
fit4.2b <- aov(RW ~ K4, data= nico_morpholin) #ancho del rostro
summary(fit4.2b)
# Df Sum Sq Mean Sq F value Pr(>F)   
#K4           3  4.106   1.369   6.548 0.0012 **
#Residuals   36  7.525   0.209 
TukeyHSD(fit4.2b)
#                     diff        lwr        upr     p adj
#North1-Central -0.9577193 -1.7226755 -0.1927631 0.0092963 **
#North2-Central -0.6250526 -1.1060967 -0.1440086 0.0066139 **
summary(fit4.2a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  2.185  0.7285   1.808  0.163
#Residuals   36 14.505  0.4029
TukeyHSD(fit4.2a)
summary(fit4.2)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K4           3  11.90   3.966   3.893 0.0165 *
#Residuals   36  36.67   1.019                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(fit4.2) 
#                     diff        lwr         upr     p adj
#North1-Central -1.4035088 -3.0921703  0.2851528 0.1322575
#North2-Central -1.1758421 -2.2377599 -0.1139243 0.0251550 *
#South-Central  -0.6418421 -1.7874310  0.5037468 0.4428400
#North2-North1   0.2276667 -1.5616193  2.0169526 0.9859303
#South-North1    0.7616667 -1.0785096  2.6018429 0.6829316
#South-North2    0.5340000 -0.7553178  1.8233178 0.6825040


fit4.3 <- aov(RL+RW ~ Clado, data= nico_morpholin) 
#            Df Sum Sq Mean Sq F value  Pr(>F)   
#Region       2  11.78   5.889   5.923 0.00587 **
#Residuals   37  36.79   0.994                   
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(fit4.3)
fit4.3a <- aov(RL ~ Clado, data= nico_morpholin)
fit4.3b <- aov(RW ~ Clado, data= nico_morpholin)
summary(fit4.3b)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#Clado        2  3.851  1.9253   9.156 0.000588 ***
#Residuals   37  7.780  0.2103 
TukeyHSD(fit4.3b)
#                    diff        lwr        upr     p adj
#North-Central -0.7018219 -1.1047878 -0.2988559 0.0003972 ***
summary(fit4.3a)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#Clado        2   2.16  1.0800    2.75  0.077 .
#Residuals   37  14.53  0.3927 
TukeyHSD(fit4.3a)
TukeyHSD(fit4.3) 
#                    diff        lwr        upr     p adj
#North-Central -1.2283806 -2.1046424 -0.3521188 0.0042629
#South-Central -0.6418421 -1.6678896  0.3842054 0.2900786
#South-North    0.5865385 -0.5074188  1.6804957 0.3992310


fit4.4 <- aov(RL+RW ~ K5, data= nico_morpholin)
fit4.4a <- aov(RL ~ K5, data= nico_morpholin)
fit4.4b <- aov(RW ~ K5, data= nico_morpholin)
summary(fit4.4b)
#            Df Sum Sq Mean Sq F value  Pr(>F)   
#K5           4  4.366  1.0915   5.258 0.00201 **
#Residuals   35  7.265  0.2076 
TukeyHSD(fit4.4b)
#North1-Center1  -0.8344444 -1.7076788  0.03878992 0.0668606 .
#North1-Center2  -1.0686667 -1.9309165 -0.20641682 0.0089550 **
#North2-Center2  -0.7360000 -1.3217834 -0.15021658 0.0078584 **
summary(fit4.4a)
#           Df Sum Sq Mean Sq F value Pr(>F)   
#K5           4  5.895  1.4738   4.779 0.0035 **
#Residuals   35 10.794  0.3084  
TukeyHSD(fit4.4a)
#                    diff        lwr         upr     p adj
#Center2-Center1  0.88500  0.1513828  1.61861720 0.0115049 *
#North2-Center2  -0.97000 -1.6840505 -0.25594945 0.0035326 **
#South1-Center2  -0.68875 -1.4461150  0.06861498 0.0894110 .
summary(fit4.4)
#           Df Sum Sq Mean Sq F value  Pr(>F)   
#K5           4  17.83   4.458   5.076 0.00247 **
#Residuals   35  30.73   0.878
TukeyHSD(fit4.4)
#                       diff        lwr        upr     p adj
#Center2-Center1  1.11922222 -0.1186749  2.3571193 0.0923550
#North1-Center1  -0.81444444 -2.6105744  0.9816855 0.6907191
#North2-Center1  -0.58677778 -1.8246749  0.6511193 0.6547130
#South1-Center1  -0.05277778 -1.3619212  1.2563656 0.9999571
#North1-Center2  -1.93366667 -3.7072029 -0.1601305 0.0268803 *
#North2-Center2  -1.70600000 -2.9108806 -0.5011194 0.0022253 **
#South1-Center2  -1.17200000 -2.4499688  0.1059688 0.0852654
#North2-North1    0.22766667 -1.5458695  2.0012029 0.9958579
#South1-North1    0.76166667 -1.0623119  2.5856452 0.7509929
#South1-North2    0.53400000 -0.7439688  1.8119688 0.7505621



fit5 <- aov(AreL+AreW ~ Clado, data= nico_morpholin)
summary(fit5)
#           Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2   1.23  0.6166   0.327  0.723
#Residuals   37  69.87  1.8883
fit5a <- aov(AreL ~ K5, data= nico_morpholin)
summary(fit5a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   1.78  0.4451   0.499  0.737
#Residuals   35  31.24  0.8925
fit5b <- aov(AreW ~ K5, data= nico_morpholin)
summary(fit5b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   0.80  0.1999   0.553  0.698
#Residuals   35  12.65  0.3614
TukeyHSD(fit5)
#                    diff       lwr       upr     p adj
#North-Central -0.2510526 -1.458645 0.9565402 0.8681277
#South-Central -0.4435526 -1.857568 0.9704627 0.7260239
#South-North   -0.1925000 -1.700103 1.3151030 0.9479140


fit5.1 <- aov(AreL+AreW ~ K4, data= nico_morpholin)
fit5.1a <- aov(AreL ~ K4, data= nico_morpholin)
fit5.1b <- aov(AreW ~ K4, data= nico_morpholin)
summary(fit5.1b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  0.798  0.2659   0.757  0.526
#Residuals   36 12.650  0.3514
summary(fit5.1a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3  1.423  0.4745   0.541  0.658
#Residuals   36 31.593  0.8776 
summary(fit5.1)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K4           3   4.12   1.373   0.738  0.537
#Residuals   36  66.98   1.861  
TukeyHSD(fit5.1)
#                       diff       lwr      upr     p adj
#North1-Central -1.111052632 -3.393399 1.171293 0.5620676
#North2-Central  0.006947368 -1.428310 1.442205 0.9999992
#South-Central  -0.443552632 -1.991897 1.104792 0.8667509
#North2-North1   1.118000000 -1.300347 3.536347 0.6029355
#South-North1    0.667500000 -1.819629 3.154629 0.8873453
#South-North2   -0.450500000 -2.193104 1.292104 0.8978448


fit5.2 <- aov(AreL+AreW ~ K5, data= nico_morpholin)
summary(fit5.2)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   4.42   1.105    0.58  0.679
#Residuals   35  66.68   1.905
fit5.2a <- aov(AreL ~ K5, data= nico_morpholin)
summary(fit5.2a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   1.78  0.4451   0.499  0.737
#Residuals   35  31.24  0.8925
fit5.2b <- aov(AreW ~ K5, data= nico_morpholin)
summary(fit5.2b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4   0.80  0.1999   0.553  0.698
#Residuals   35  12.65  0.3614 
TukeyHSD(fit5.2)
#                      diff       lwr      upr     p adj
#Center2-Center1 -0.2534444 -2.076776 1.569887 0.9943736
#North1-Center1  -1.2444444 -3.890012 1.401123 0.6611111
#North2-Center1  -0.1264444 -1.949776 1.696887 0.9996292
#South1-Center1  -0.5769444 -2.505217 1.351328 0.9093440
#North1-Center2  -0.9910000 -3.603289 1.621289 0.8100835
#North2-Center2   0.1270000 -1.647701 1.901701 0.9995802
#South1-Center2  -0.3235000 -2.205855 1.558855 0.9874046
#North2-North1    1.1180000 -1.494289 3.730289 0.7339304
#South1-North1    0.6675000 -2.019087 3.354087 0.9518151
#South1-North2   -0.4505000 -2.332855 1.431855 0.9577575


fit5.3 <- aov(AreL+AreW ~ COI, data= nico_morpholin)
summary(fit5.3)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5   9.93   1.987   1.104  0.376
#Residuals   34  61.17   1.799
fit5.3a <- aov(AreL ~ COI, data= nico_morpholin)
summary(fit5.3a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5  2.905  0.5811   0.656  0.659
#Residuals   34 30.111  0.8856
fit5.3b <- aov(AreW ~ COI, data= nico_morpholin)
summary(fit5.3b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5  2.456  0.4912   1.519   0.21
#Residuals   34 10.992  0.3233 
TukeyHSD(fit5.3)
#                      diff       lwr      upr     p adj
#Center2-Center1 -0.2534444 -2.113537 1.606648 0.9983538
#North1-Center1  -1.2444444 -3.943350 1.454461 0.7315585
#North2-Center1  -0.1264444 -1.986537 1.733648 0.9999449
#South1-Center1  -1.4069444 -3.839705 1.025816 0.5127154
#South2-Center1   0.2530556 -2.179705 2.685816 0.9995538
#North1-Center2  -0.9910000 -3.655956 1.673956 0.8686786
#North2-Center2   0.1270000 -1.683481 1.937481 0.9999356
#South1-Center2  -1.1535000 -3.548541 1.241541 0.6948256
#South2-Center2   0.5065000 -1.888541 2.901541 0.9871978
#North2-North1    1.1180000 -1.546956 3.782956 0.8007921
#South1-North1   -0.1625000 -3.254485 2.929485 0.9999846
#South2-North1    1.4975000 -1.594485 4.589485 0.6899049
#South1-North2   -1.2805000 -3.675541 1.114541 0.5955238
#South2-North2    0.3795000 -2.015541 2.774541 0.9966262
#South2-South1    1.6600000 -1.202622 4.522622 0.5098001

fit5.4 <- aov(AreL+AreW ~ BFD_gdi, data= nico_morpholin)
summary(fit5.4)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6   9.98   1.664   0.898  0.508
#Residuals   33  61.12   1.852 
fit5.4a <- aov(AreL ~ BFD_gdi, data= nico_morpholin)
summary(fit5.4a)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  3.021  0.5035   0.554  0.763
#Residuals   33 29.995  0.9089
fit5.4b <- aov(AreW ~ BFD_gdi, data= nico_morpholin)
summary(fit5.4b)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6   2.47  0.4116   1.237  0.313
#Residuals   33  10.98  0.3327 
TukeyHSD(fit5.4)
#                      diff       lwr      upr     p adj
#Center2-Center1 -0.2534444 -2.215048 1.708159 0.9996024
#North1-Center1  -1.2444444 -4.090638 1.601749 0.8124873
#North2-Center1  -0.1264444 -2.088048 1.835159 0.9999933
#South1-Center1  -1.4069444 -3.972468 1.158580 0.6078946
#South2-Center1   0.4455556 -4.054671 4.945782 0.9999153
#South3-Center1   0.1888889 -2.657304 3.035082 0.9999920
#North1-Center2  -0.9910000 -3.801391 1.819391 0.9217164
#North2-Center2   0.1270000 -1.782285 2.036285 0.9999919
#South1-Center2  -1.1535000 -3.679246 1.372246 0.7804562
#South2-Center2   0.6990000 -3.778669 5.176669 0.9988304
#South3-Center2   0.4423333 -2.368057 3.252724 0.9987750
#North2-North1    1.1180000 -1.692391 3.928391 0.8698828
#South1-North1   -0.1625000 -3.423224 3.098224 0.9999985
#South2-North1    1.6900000 -3.239751 6.619751 0.9309065
#South3-North1    1.4333333 -2.052527 4.919194 0.8516152
#South1-North2   -1.2805000 -3.806246 1.245246 0.6890436
#South2-North2    0.5720000 -3.905669 5.049669 0.9996275
#South3-North2    0.3153333 -2.495057 3.125724 0.9998242
#South2-South1    1.8525000 -2.920711 6.625711 0.8823192
#South3-South1    1.5958333 -1.664891 4.856557 0.7221632
#South3-South2   -0.2566667 -5.186418 4.673085 0.9999981



fit6 <- aov(AW ~ Clado, data= nico_morpholin)
summary(fit6)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Clado        2   14.4   7.198   2.159   0.13
#Residuals   37  123.4   3.334 
TukeyHSD(fit6)
#                    diff       lwr       upr     p adj
#North-Central -0.8943725 -2.499025 0.7102800 0.3715547
#South-Central -1.4967763 -3.375723 0.3821708 0.1405010
#South-North   -0.6024038 -2.605711 1.4009029 0.7449325

fit6.1 <- aov(AW ~ K4, data= nico_morpholin)
summary(fit6.1)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#K4           3  25.27   8.422   2.695 0.0604 .
#Residuals   36 112.50   3.125
TukeyHSD(fit6.1)
#                     diff        lwr       upr     p adj
#North1-Central -2.5638596 -5.5216540 0.3939347 0.1090243
#North2-Central -0.3935263 -2.2535405 1.4664879 0.9404159
#South-Central  -1.4967763 -3.5033455 0.5097929 0.2036241
#North2-North1   2.1703333 -0.9637108 5.3043775 0.2608571
#South-North1    1.0670833 -2.1560983 4.2902650 0.8091968
#South-North2   -1.1032500 -3.3615691 1.1550691 0.5591991

fit6.2 <- aov(AW ~ K5, data= nico_morpholin)
summary(fit6.2)
#            Df Sum Sq Mean Sq F value Pr(>F)
#K5           4  25.69   6.422   2.006  0.115
#Residuals   35 112.08   3.202 
TukeyHSD(fit6.2)
#                      diff       lwr       upr     p adj
#Center2-Center1 -0.2987778 -2.662645 2.0650894 0.9960987
#North1-Center1  -2.7211111 -6.150970 0.7087479 0.1752249
#North2-Center1  -0.5507778 -2.914645 1.8130894 0.9615861
#South1-Center1  -1.6540278 -4.153946 0.8458901 0.3350497
#North1-Center2  -2.4223333 -5.809048 0.9643811 0.2615526
#North2-Center2  -0.2520000 -2.552819 2.0488194 0.9977628
#South1-Center2  -1.3552500 -3.795637 1.0851375 0.5093777
#North2-North1    2.1703333 -1.216381 5.5570478 0.3664908
#South1-North1    1.0670833 -2.415955 4.5501217 0.9020386
#South1-North2   -1.1032500 -3.543637 1.3371375 0.6930770

fit6.3 <- aov(AW ~ COI, data= nico_morpholin)
summary(fit6.3)
#            Df Sum Sq Mean Sq F value Pr(>F)
#COI          5  29.91   5.982   1.886  0.123
#Residuals   34 107.86   3.172
TukeyHSD(fit6.3)
#                      diff       lwr       upr     p adj
#Center2-Center1 -0.2987778 -2.768747 2.1711912 0.9990715
#North1-Center1  -2.7211111 -6.304919 0.8626966 0.2252633
#North2-Center1  -0.5507778 -3.020747 1.9191912 0.9837678
#South1-Center1  -2.3802778 -5.610678 0.8501229 0.2535708
#South2-Center1  -0.9277778 -4.158178 2.3026229 0.9518260
#North1-Center2  -2.4223333 -5.961060 1.1163933 0.3284652
#North2-Center2  -0.2520000 -2.656091 2.1520913 0.9995368
#South1-Center2  -2.0815000 -5.261814 1.0988139 0.3767482
#South2-Center2  -0.6290000 -3.809314 2.5513139 0.9905511
#North2-North1    2.1703333 -1.368393 5.7090599 0.4483046
#South1-North1    0.3408333 -3.764934 4.4466009 0.9998523
#South2-North1    1.7933333 -2.312434 5.8991009 0.7730111
#South1-North2   -1.8295000 -5.009814 1.3508139 0.5184883
#South2-North2   -0.3770000 -3.557314 2.8033139 0.9991581
#South2-South1    1.4525000 -2.348702 5.2537021 0.8552674

fit6.4 <- aov(AW ~ BFD_gdi, data= nico_morpholin)
summary(fit6.4)
#            Df Sum Sq Mean Sq F value Pr(>F)
#BFD_gdi      6  30.07   5.012   1.536  0.197
#Residuals   33 107.69   3.263
TukeyHSD(fit6.4)
#                      diff       lwr      upr     p adj
#Center2-Center1 -0.2987778 -2.902619 2.305064 0.9997998
#North1-Center1  -2.7211111 -6.499162 1.056939 0.2932248
#North2-Center1  -0.5507778 -3.154619 2.053064 0.9937457
#South1-Center1  -2.3802778 -5.785767 1.025211 0.3264928
#South2-Center1  -0.5777778 -6.551400 5.395845 0.9999261
#South3-Center1  -1.0444444 -4.822495 2.733606 0.9750001
#North1-Center2  -2.4223333 -6.152859 1.308193 0.4119464
#North2-Center2  -0.2520000 -2.786393 2.282393 0.9999132
#South1-Center2  -2.0815000 -5.434187 1.271187 0.4651698
#South2-Center2  -0.2790000 -6.222679 5.664679 0.9999990
#South3-Center2  -0.7456667 -4.476193 2.984859 0.9953927
#North2-North1    2.1703333 -1.560193 5.900859 0.5414593
#South1-North1    0.3408333 -3.987467 4.669134 0.9999780
#South2-North1    2.1433333 -4.400442 8.687109 0.9437899
#South3-North1    1.6766667 -2.950481 6.303815 0.9118664
#South1-North2   -1.8295000 -5.182187 1.523187 0.6133064
#South2-North2   -0.0270000 -5.970679 5.916679 1.0000000
#South3-North2   -0.4936667 -4.224193 3.236859 0.9995438
#South2-South1    1.8025000 -4.533483 8.138483 0.9712091
#South3-South1    1.3358333 -2.992467 5.664134 0.9573849
#South3-South2   -0.4666667 -7.010442 6.077109 0.9999878

####################################################


S <- cov(nico_morpholin[,12:25])
S

sum(diag(S))
s.eigen <- eigen(S)
s.eigen

for (s in s.eigen$values) {
  print(s / sum(s.eigen$values))
}
#The first two principal components account for 94% of the total variance
#A scree graph of the eigenvalues can be plotted to visualize the proportion of variance explained by each subsequential eigenvalue.

plot(s.eigen$values, xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph')
lines(s.eigen$values)
s.eigen$vectors

#Principal Component Analysis with R
#Computing the principal components in R is straightforward with the functions 
#prcomp() and princomp(). The difference between the two is simply the method 
#employed to calculate PCA. According to prcomp:
#The calculation is done by a singular value decomposition of the (centered and
# possibly scaled) data matrix, not by using eigen on the covariance matrix. 
# This is generally the preferred method for numerical accuracy.

nico.pca <-  prcomp(nico_morpholin[,12:25])
nico.pca

#Although we didn't use the preferred method of applying singular value 
#decomposition, the components reported by the prcomp() are the same as what was computed earlier save arbitrary scalings of ???1 to some of the eigenvectors.
#The summary method of prcomp() also outputs the proportion of variance 
#explained by the components.

summary(nico.pca)

#Plotting of Principal Components
#The first two principal components are often plotted as a scatterplot which may reveal interesting features of the data, such as departures from normality, outliers or non-linearity. The first two principal components are evaluatedfor each observation vector and plotted.

#The ggfortify package provides a handy method for plotting the first two 
#principal components with autoplot().

library(ggfortify)

pca.plot <- autoplot(nico.pca, data = nico_morpholin, colour = 'BFD_gdi', frame = TRUE)
pca.plot
summary(pca.plot)

pca.plot_1 <- autoplot(nico.pca, data = nico_morpholin, colour = 'COI', frame = TRUE)
pca.plot_1

pca.plot_2 <- autoplot(nico.pca, data = nico_morpholin, colour = 'K4', frame = TRUE)
pca.plot_2

pca.plot_3 <- autoplot(nico.pca, data = nico_morpholin, colour = 'K5', frame = TRUE)
pca.plot_3

pca.plot_4 <- autoplot(nico.pca, data = nico_morpholin, colour = 'Clado', frame = TRUE)
pca.plot_4

ggarrange(pca.plot_4, pca.plot_2, pca.plot_3, pca.plot_1, pca.plot, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)

#Interpreting Principal Components
#Interpretation of principal components is still a heavily researched topic in statistics, and although the components may be readily interpreted in most settings, this is not always the case (Joliffe, 2002).
#
#One method of interpretation of the principal components is to calculate the correlation between the original data and the component. The autoplot() function also generates a nice data table with the original variables and the calculated PCs, which we will use here to find the correlations.
#First, compute the correlations between the data and the calculated components of the covariance matrix S.


##S2 analysis PCA
##
S2 <-cov(nico_morpholin[,23:25])
S2

sum(diag(S2))
s.eigen <- eigen(S2)
s.eigen


for (s in s.eigen$values) {
  print(s / sum(s.eigen$values))
}
#The first two principal components account for 94% of the total variance
#A scree graph of the eigenvalues can be plotted to visualize the proportion of variance explained by each subsequential eigenvalue.

plot(s.eigen$values, xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph')
lines(s.eigen$values)
s.eigen$vectors

S2pca <-  prcomp(nico_morpholin[,23:25])
class(S2pca)
summary(S2pca)
S2pca$x
class(pcaS2_2)

library(ggplot2)
library(ggfortify)
library(gridExtra)
library(ggpubr)

pcaS2_1f <- autoplot(S2pca, data = nico_morpholin, colour = 'Clado', frame = TRUE)
pcaS2_1f

pcaS2_2f <- autoplot(S2pca, data = nico_morpholin, colour = 'K4', frame = TRUE)
pcaS2_2 <- autoplot(S2pca, data = nico_morpholin, colour = 'K4')
pcaS2_2f

pcaS2_3f <- autoplot(S2pca, data = nico_morpholin, colour = 'K5', frame = TRUE)
pcaS2_3 <- autoplot(S2pca, data = nico_morpholin, colour = 'K5')
pcaS2_3f

pcaS2_4f <- autoplot(S2pca, data = nico_morpholin, colour = 'COI', frame = TRUE)
pcaS2_4 <- autoplot(S2pca, data = nico_morpholin, colour = 'COI')
pcaS2_4f

pcaS2_5f <- autoplot(S2pca, data = nico_morpholin, colour = 'BFD_gdi', frame = TRUE)
pcaS2_5 <- autoplot(S2pca, data = nico_morpholin, colour = 'BFD_gdi')
pcaS2_5f

ggarrange(pcaS2_1f, pcaS2_2f, pcaS2_3f, pcaS2_4f, pcaS2_5f, 
           labels = c("A", "B", "C", "D", "E"),
           ncol = 3, nrow = 2)

#################
####Rostro
rostro <-cov(nico_morpholin[,18:19])
rostro

sum(diag(rostro))
s.eigen <- eigen(rostro)
s.eigen


for (s in s.eigen$values) {
  print(s / sum(s.eigen$values))
}


Rpca <-  prcomp(nico_morpholin[,18:19])
Rpca
summary(Rpca)

pcaR_1 <- autoplot(Rpca, data = nico_morpholin, colour = 'BFD_gdi', frame = TRUE)
pcaR_1

pcaR_2 <- autoplot(Rpca, data = nico_morpholin, colour = 'COI', frame = TRUE)
pcaR_2

pcaR_3 <- autoplot(Rpca, data = nico_morpholin, colour = 'K4', frame = TRUE)
pcaR_3

pcaR_4 <- autoplot(Rpca, data = nico_morpholin, colour = 'K5', frame = TRUE)
pcaR_4

pcaR_5 <- autoplot(Rpca, data = nico_morpholin, colour = 'Clado', frame = TRUE)
pcaR_5

ggarrange(pcaR_5, pcaR_3, pcaR_4, pcaR_2, pcaR_1, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)


##################
#quela
########

Qpca <-  prcomp(nico_morpholin[,20:22])
Qpca
summary(Qpca)

pcaQ_1 <- autoplot(Qpca, data = nico_morpholin, colour = 'Clado', frame = TRUE)
pcaQ_1

pcaQ_2 <- autoplot(Qpca, data = nico_morpholin, colour = 'K4', frame = TRUE)
pcaQ_2

pcaQ_3 <- autoplot(Qpca, data = nico_morpholin, colour = 'K5', frame = TRUE)
pcaQ_3

pcaQ_4 <- autoplot(Qpca, data = nico_morpholin, colour = 'COI', frame = TRUE)
pcaQ_4

pcaQ_5 <- autoplot(Qpca, data = nico_morpholin, colour = 'BFD_gdi', frame = TRUE)
pcaQ_5

ggarrange(pcaQ_1, pcaQ_2, pcaQ_3, pcaQ_4, pcaQ_5, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)


#Total length

TLpca <-  prcomp(nico_morpholin[,12:14])
TLpca
summary(TLpca)

pcaL <- autoplot(TLpca, data = nico_morpholin, colour = 'Clado', frame = TRUE)
pcaL

pcaL2 <- autoplot(TLpca, data = nico_morpholin, colour = 'K4', frame = TRUE)
pcaL2

pcaL3 <- autoplot(TLpca, data = nico_morpholin, colour = 'K5', frame = TRUE)
pcaL3

pcaL4 <- autoplot(TLpca, data = nico_morpholin, colour = 'COI', frame = TRUE)
pcaL4

pcaL5 <- autoplot(TLpca, data = nico_morpholin, colour = 'BFD_gdi', frame = TRUE)
pcaL5

ggarrange(pcaL, pcaL2, pcaL3, pcaL4, pcaL5, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)

#################
####Areola
areola <-cov(nico_morpholin[,16:17])
areola

sum(diag(areola))
s.eigen <- eigen(areola)
s.eigen


for (s in s.eigen$values) {
  print(s / sum(s.eigen$values))
}


Apca <-  prcomp(nico_morpholin[,16:17])
Apca
summary(Apca)

pcaA_1 <- autoplot(Apca, data = nico_morpholin, colour = 'Clado', frame = TRUE)
pcaA_1

pcaA_2 <- autoplot(Apca, data = nico_morpholin, colour = 'K4', frame = TRUE)
pcaA_2

pcaA_3 <- autoplot(Apca, data = nico_morpholin, colour = 'K5', frame = TRUE)
pcaA_3

pcaA_4 <- autoplot(Apca, data = nico_morpholin, colour = 'COI', frame = TRUE)
pcaA_4

pcaA_5 <- autoplot(Apca, data = nico_morpholin, colour = 'BFD_gdi', frame = TRUE)
pcaA_5

ggarrange(pcaA_1, pcaA_2, pcaA_3, pcaA_4, pcaA_5, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)

#################
####abdomen ancho
AW <-cov(nico_morpholin[15])
AW

sum(diag(AW))
s.eigen <- eigen(AW)
s.eigen


for (s in s.eigen$values) {
  print(s / sum(s.eigen$values))
}

AWpca <-  prcomp(nico_morpholin$AW)
AWpca
summary(AWpca)

pcaAW_1 <- autoplot(AWpca, data = nico_morpholin, colour = 'Clado', frame = TRUE)
pcaAW_1



pcaA_2 <- autoplot(Apca, data = nico_morpholin, colour = 'K4', frame = TRUE)
pcaA_2

pcaA_3 <- autoplot(Apca, data = nico_morpholin, colour = 'K5', frame = TRUE)
pcaA_3

pcaA_4 <- autoplot(Apca, data = nico_morpholin, colour = 'COI', frame = TRUE)
pcaA_4

pcaA_5 <- autoplot(Apca, data = nico_morpholin, colour = 'BFD_gdi', frame = TRUE)
pcaA_5

ggarrange(pcaA_1, pcaA_2, pcaA_3, pcaA_4, pcaA_5, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)


citation()
