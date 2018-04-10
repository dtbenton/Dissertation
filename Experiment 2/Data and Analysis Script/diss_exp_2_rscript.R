########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 2 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################
# load all relevant libraries:
library(lme4)
library(nlme)
library(boot)
library(car) 
library(reshape2)
library(ggplot2)
library(ez)
library(plyr)
library(ggsignif)
library(rcompanion)
options(scipen=9999)

# load data
D = read.csv(file.choose(), header = TRUE)
D = as.data.frame(D[1:64,])

# reorder columns
D = as.data.frame(D[,c(1,2,3,7,5,9,4,8,6,10,11,12)])

# reshape the data
D_tall = reshape(D, varying = 3:10, v.names = "measure", 
                 timevar = "condition", idvar = "ID", 
                 new.row.names = 1:512, direction = "long")

# order data
D_tall = D_tall[order(D_tall$ID),]