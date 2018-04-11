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
options(scipen=9999)

# load data
D = read.csv(file.choose(), header = TRUE)

# reorder columns
D = as.data.frame(D[,c(1,2,3,7,5,9,4,8,6,10,11,12)])

# reshape the data
D_tall = reshape(D, varying = 3:10, v.names = "measure", 
                 timevar = "condition", idvar = "ID", 
                 new.row.names = 1:512, direction = "long")

# order data
D_tall = D_tall[order(D_tall$ID),]


# add q.type.cat column
D_tall$q.type.cat = as.factor(rep(c(1:2), each = 4, times = 64))


# set appropriate factor variables in "tall" data
D_tall$condition = as.factor(D_tall$condition)
D_tall$group = as.factor(D_tall$group)
D_tall$q.type = as.factor(D_tall$q.type)
D_tall$q.type.cat = as.factor(D_tall$q.type.cat)
D_tall$measure.2 = (100-D_tall$measure)



########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################


## NORMALITY CHECKS

par(mfrow=c(4,2)) 
for (ii in 1:8)  hist(D_tall$measure[D_tall$condition==ii], breaks=5)
par(mfrow=c(1,1)) 

# formal test of normality
shapiro.ps = rep(0,8)
for(i in 1:8) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$condition==i])
  shapiro.ps[i] = shap.calc$p.value
}


## HOMOSKEDASTICITY CHECKS
# plot the boxplots
# perceptual
boxplot(D_tall$measure[D_tall$condition==c(1:4)]~D_tall$condition[D_tall$condition==c(1:4)])

# causal
boxplot(D_tall$measure[D_tall$condition==c(5:8)]~D_tall$condition[D_tall$condition==c(5:8)])

# formal test of equal variance
# perceptual
leveneTest(D_tall$measure[D_tall$condition==c(1:4)], as.factor(D_tall$condition[D_tall$condition==c(1:4)]), center=median) # used 'median' because it's a better measure of central tendency given the non-normality

# causal 
leveneTest(D_tall$measure[D_tall$condition==c(5:8)], as.factor(D_tall$condition[D_tall$condition==c(5:8)]), center=median) # used 'median' because it's a better measure of central tendency given the non-normality


## ASSUMPTION CHECK NOTES ##
# Despite the fact that there is no evidence of heteroskedasticity, because
# there is evidence of non-normality, non-parametric boostrapping and permuation
# testing will be used to estimate confidence intervals and for hypothesis checking


########################################################
########################################################
########################################################
#############                              #############
#############            Models            #############
#############                              #############
########################################################
########################################################
########################################################



##############################
#### PRELIMINARY ANALYSIS ####
##############################
# analysis to determine effect of question type and 
lme.fit.prelim = lme(measure~q.type.cat+group+q.type.cat:group, random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)

## post hocs ##

# red circle first
set.seed(2018)
b = rep(0,4000) 
for(i in 1:4000){
  y = sample(D_tall$measure.2, replace=TRUE)
  lm_1 = lmer(y[D_tall$group==0]~D_tall$q.type.cat[D_tall$group==0]+(1|D_tall$ID[D_tall$group==0]), 
              data=D_tall)
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lmer(D_tall$measure.2[D_tall$group==0]~D_tall$q.type.cat[D_tall$group==0]+(1|D_tall$ID[D_tall$group==0]), 
              data=D_tall)
beta_actual = fixed.effects(lm.fit)[[2]]
c(beta_actual, sum(abs(b) > beta_actual)/4000, sum(abs(b) < beta_actual)/4000,
  sum(b > beta_actual)/4000, sum(b < beta_actual)/4000)



# blue circle first
set.seed(2018)
b = rep(0,4000) 
for(i in 1:4000){
  y = sample(D_tall$measure.2, replace=TRUE)
  lm_1 = lmer(y[D_tall$group==1]~D_tall$q.type.cat[D_tall$group==1]+(1|D_tall$ID[D_tall$group==1]), 
              data=D_tall)
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lmer(D_tall$measure.2[D_tall$group==1]~D_tall$q.type.cat[D_tall$group==1]+(1|D_tall$ID[D_tall$group==1]), 
              data=D_tall)
beta_actual = fixed.effects(lm.fit)[[2]]
c(beta_actual, sum(abs(b) > beta_actual)/4000, sum(abs(b) < beta_actual)/4000,
  sum(b > beta_actual)/4000, sum(b < beta_actual)/4000)





## PRELIMINARY ANALYSIS NOTES
# no effect of question type (perceptual first vs causal first) or location (red first or blue 
# first in training sequence)
condition_barplot = ggplot(D_tall, aes(group, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~group) + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() 
