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





#######################
#### MAIN ANALYSIS ####
#######################
# perceptual question
sub.percep = subset(D_tall, ! condition %in% c(5:8))
lme.fit.main.percep = lme(measure.2~as.factor(condition), random=~1|ID, data = sub.percep)
anova.lme(lme.fit.main.percep)



# causal question
sub.causal = subset(D_tall, ! condition %in% c(1:4))
lme.fit.main.causal = lme(measure.2~condition, random=~1|ID, data = sub.causal)
anova.lme(lme.fit.main.causal)




###################################
## FOLLOW-UP PLANNED COMPARISONS ##
###################################
#  permutation global function
perm_func = function(p,v){
  set.seed(2018)
  b = rep(0,4000) 
  for(i in 1:4000){
    x = factor(D_tall$condition, levels=c(p,v)) 
    y = sample(D_tall$measure.2, replace=TRUE)
    lm_1 = lmer(y ~ x + (1|ID), data=D_tall) 
    b[i] = fixed.effects(lm_1)[2]
  }
  
  lm.fit = lmer(D_tall$measure.2~factor(D_tall$condition, levels=c(p,v))+(1|ID), data=D_tall)
  beta_actual = fixed.effects(lm.fit)[[2]]
  c(beta_actual, sum(abs(b) > beta_actual)/4000, sum(abs(b) < beta_actual)/4000,
    sum(b > beta_actual)/4000, sum(b < beta_actual)/4000)
}



# bootstrap global function
boot_mean = as.data.frame(matrix(NA, nrow=8, ncol=3, byrow=TRUE))
for(i in 1:nrow(boot_mean)){ # want number of iterations to equal number of rows, especially because we're filling in by row
  set.seed(2018)
  boot_func = function(data,b,formula, p){ 
    d= data[b,] 
    x = d[d$condition==i,8]
    dif.1 =  mean(x, data=D_tall) 
    return(dif.1)
  }
  
  GBGR_P_boot = boot(D_tall, boot_func, R=4000) 
  boot_mean[i,] = c(GBGR_P_boot$t0, GBGR_P_boot$t0  + 1.96*-sd(GBGR_P_boot$t), 
                    GBGR_P_boot$t0  + 1.96*sd(GBGR_P_boot$t))
}



############################
## PERM TESTS: PERCEPTUAL ##
############################
# GBGR-P V. GBgapGR-P
perm_func(1,2)

# GBGR-P V. GBRG-P
perm_func(1,3)

# GBGR-P V. GBgapRG-P
perm_func(1,4)

# GBgapGR-P V. GBRG-P
perm_func(2,3)

# GBRG-P V. GBgapRG-P
perm_func(2,4)



########################
## PERM TESTS: CAUSAL ##
########################
# GBGR-C V. GBRG-C
perm_func(5,6)

# GBGR-C V. GBRG-C
perm_func(5,7)

# GBGR-C V. GBgapRG-C
perm_func(5,8)

# GBgapGR-C V. GBRG-C
perm_func(6,7)

# GBRG-C V. GBgapRG-C
perm_func(7,8)


# GBgapGR vs GBgapRG
perm_func(6,8)



#####################################
## BOOT TESTS: PERCEPTUAL & CAUSAL ##
#####################################
V1       V2       V3
1 22.34375 16.33642 28.35108 # P GBGR
2 49.54688 43.74369 55.35006 # P GBgapGR
3 22.03125 16.32894 27.73356 # P GBRG
4 47.62500 42.26285 52.98715 # P GBgapRG
5 17.37500 11.92218 22.82782 # C GBGR
6 49.03125 43.31298 54.74952 # C GBgapGR
7 15.34375 10.72973 19.95777 # C GBRG
8 42.84375 36.81343 48.87407 # C GBgapRG





########################################################
########################################################
########################################################
#############                              #############
#############            Figures           #############
#############                              #############
########################################################
########################################################
########################################################

# Create 'F_tall' data frame to use for ggplot
F_tall = D_tall

# rename levels of 'condition' and 'q.type.cat' factors
F_tall$condition = revalue(x = as.factor(F_tall$condition), 
                           c("1" = "GBGR-P", "2"="GBgapGR-P", "3" = "GBRG-P", 
                             "4" = "GBgapRG-P", "5" = "GBGR-C", "6" = "GBgapGR-C",
                             "7" = "GBRG-C", "8" = "GBgapRG-C"))
F_tall$q.type.cat = revalue(x = as.factor(F_tall$q.type.cat), 
                            c("1" = "Perceptual Question", "2"="Causal Question"))

F_tall$condition.2 = rep(c(1:4), times = 128)
F_tall$condition.2 = revalue(x = as.factor(F_tall$condition.2), 
                           c("1" = "GBGR", "2"="GBgapGR", "3" = "GBRG", 
                             "4" = "GBgapRG"))

F_tall$q.type.cat.2 = rep(c(0,1), each = 4, times = 64)
F_tall$q.type.cat.2 = revalue(x = as.factor(F_tall$q.type.cat.2), 
                              c("0" = "Perceptual Question", "1"="Causal Question"))


# OMNIBUS ANALYSIS FIGURE
condition_barplot = ggplot(F_tall, aes(q.type.cat.2, measure.2, fill = condition.2)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  geom_signif(annotations = c("p < .0001","p < .0001","p < .0001","p < .0001"),
              y_position = c(68,64,64,60), xmin=c(.6,1.125,.6,.9), 
              xmax=c(1.35,1.35,.9,1.1), 
              tip_length = 0.00375) +
  geom_signif(annotations = c("p < .0001","p < .0001","p < .0001","p < .0001"),
              y_position = c(68,64,64,60), xmin=c(2.275,2.275,1.875,1.875), 
              xmax=c(1.6,2,1.6,2.075), 
              tip_length = 0.00375) +
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 75)) +
  theme_classic() +
  scale_fill_manual(values = c("white", "gray81", "gray38", "black")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials") # change the main x-axis label

#########################################
#### INDIVIDUAL DIFFERENCES ANALYSIS ####
#########################################
# INDIVIDUAL DIFFERENCE PLOTS FOR BOTH THE PERCEPTUAL AND CAUSAL QUESTIONS
condition_barplot = ggplot(F_tall, aes(q.type.cat.2, measure.2, fill = condition.2)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  scale_fill_manual(values = c("white", "gray81", "gray38", "black")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x = "Test trials")

# INDIVIDUAL DIFFERENCE PLOTS FOR PERCEPTUAL QUESTION
condition_barplot = ggplot(sub.percep, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  scale_fill_manual(values=c("#000000", "#999999")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")


# INDIVIDUAL DIFFERENCE PLOTS FOR CAUSAL QUESTION
condition_barplot = ggplot(sub.causal, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 95)) +
  scale_fill_manual(values=c("#000000", "#999999")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")




## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES ##
# install and load post-hoc chi-square test package (called follow up binomial tests)


# create contigency table with counts
chi.data = matrix(c(0,1,0,1,2,1,47,48,2,2,0,1,13,11,51,53),2)
dimnames(chi.data) = list(c("Perceptual", "Causal"), c("Markov", "Independence", "faux-Independence",
                                                       "Spatial Gap","All Familiar","All Inconsistent",
                                                       "Other","All Associative"))
## COUNTS
# Percep: M=0, I=0, fI=6, SG=47, AF=2, AI=0, O=9, AA=55
# Causal: M=1, I=1, fI=5, SG=48, AF=2, AI=1, O=7, AA=57

# the chi data 
              Markov Independence faux-Independence Spatial Gap All Familiar All Inconsistent Other All Associative
Perceptual      0            0                 2          47            2                0    13              51
Causal          1            1                 1          48            2                1    11              53



# run chisq.test() on chi table
chi.test = chisq.test(chi.data, simulate.p.value = TRUE) 


## run post-hoc tests (where you test two particular cell means): binomial tests ##
# perceptual question

# Define a function to run post-hoc binomial tests
binom_func = function(x,y){
  bin_test = binom.test(x, x+y, p = 0.5, alternative = "two.sided")
  return(bin_test$p.value)
}


                  ##
## PERCEPTUAL CONDITION POST-HOC TESTS: ##
                  ##
  
# Independence comparisons
# M v f-I
binom_func(0,2)
# I v f-I
binom_func(0,2)
# I v T
binom_func(23,4)
# I v O
binom_fun(18,12)
# I v F
binom_func(23,7)
  
  
# Markov comparisons
# M v T
binom_func(18,4)
# M v O
binom_func(18,12)
# M v IC
binom_func(18,0)
# M v AA
binom_func(18,34)  
  
  
  
  
  
  
                      ##
## PERCEPTUAL CONDITION POST-HOC TESTS: ##
                      ##

# Independence comparisons
# M v I
binom_func(18,23)
# I v T
binom_func(23,4)
# I v O
binom_fun(18,12)
# I v F
binom_func(23,7)


# Markov comparisons
# M v T
binom_func(18,4)
# M v O
binom_func(18,12)
# M v IC
binom_func(18,0)
# M v AA
binom_func(18,34)   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # LEGEND: 
# M = Markov; I = Independence; 
# faux-Independence;
# SG = Spatial Gap; 
# AF = ALL Familiar, AI = All Inconsistent
# AA = All Associative (sum of all associative options)



# check structure of data
str(chi.data)

