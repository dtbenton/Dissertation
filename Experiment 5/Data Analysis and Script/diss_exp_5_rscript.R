########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 5 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################
# load libraries:
library(lme4)
library(nlme)
library(boot)
library(car) 
library(reshape2)
library(ggplot2)
library(ez)
library(plyr)
library(ggsignif)
options(scipen=9999)

# load data
D = read.csv(file.choose(), header = TRUE)


D_tall = reshape(D, varying = c(16:19), v.names = "measure", 
                 timevar = "group", idvar = "ID", 
                 direction = "long")

D_tall = D_tall[order(D_tall$ID),]


# add median split age column
D$med.split.age = rep(0,25)
for(i in 1:25){
  compute = ifelse(D$age[i]<=18.06,0,1)
  D$med.split.age[i] = compute
}

# reorder columns
D = as.data.frame(D[,c(1:6,22,7:21)])

# set appropriate factor variables in "tall" data
D_tall$sex = as.factor(D_tall$sex)
D_tall$group = as.factor(D_tall$group)
D_tall$med.split.age = as.factor(D_tall$med.split.age)
D_tall$hab.stim.order = as.factor(D_tall$hab.stim.order)
D_tall$test.stim.order = as.factor(D_tall$test.stim.order)



########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################

## NORMALITY CHECKS

par(mfrow=c(2,2)) 
for (ii in 1:4)  hist(D_tall$measure[D_tall$group==ii], breaks=5)
par(mfrow=c(1,1)) 

# formal test of normality
shapiro.ps = rep(0,4)
for(i in 1:4) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$group==i])
  shapiro.ps[i] = shap.calc$p.value
}

