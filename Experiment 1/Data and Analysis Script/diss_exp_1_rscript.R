########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 1 SCRIPT     #############
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

# reorder columns
D = as.data.frame(D[,c(1,2,3,9,5,7,4,10,6,8,11,12)])

# reshape the data
D_tall = reshape(D, varying = 3:10, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:512, direction = "long")

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


# names of column in 'D_tall'
"ID"        "exp"       "q.type"    "group"     "condition" "measure"   "conditon" 


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
. = rep(0,8)
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
lme.fit.prelim = lme(measure~q.type+group+q.type:group, random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)

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


# omnibus analysis
lme.fit.main.omnibus = lme(measure.2~as.factor(condition)+as.factor(q.type)+as.factor(group)+as.factor(condition):as.factor(q.type), 
                           random=~1|ID, data=D_tall)
anova.lme(lme.fit.main.omnibus)



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
# GBGR-P V. GRBG-P
perm_func(1,2)

# GBGR-P V. GBRG-P
perm_func(1,3)

# GBGR-P V. GRBG-P
perm_func(1,4)

# GRGB-P V. GBRG-P
perm_func(2,3)

# GRGB-P V. GRBG-P
perm_func(2,4)



########################
## PERM TESTS: CAUSAL ##
########################
# GBGR-P V. GRBG-P
perm_func(5,6)

# GBGR-P V. GBRG-P
perm_func(5,7)

# GBGR-P V. GRBG-P
perm_func(5,8)

# GRGB-P V. GBRG-P
perm_func(6,7)

# GRGB-P V. GRBG-P
perm_func(6,8)



################
# RELEVANT BFs #
################
# define the null and alternative models #
lm.null = lme(measure.2~1, random=~1|ID, data = sub.causal)
lm.alt = lme(measure.2~as.factor(condition), random=~1|ID, data = sub.causal)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01

# transform BF01 into posterior probabilities - the percentage (i.e., posterior likelihood) of evidence in favor the null
BF01.posterior = BF01/(1+BF01)

# transform BF10 into posterior probabilities - the percentage (i.e., posterior likelihood) of evidence in favor the alternative (e.g., 4% likely under the null than the alternative)
BF01.posterior = BF10/(1+BF10)


#####################################
## BOOT TESTS: PERCEPTUAL & CAUSAL ##
#####################################
V1       V2       V3
1 26.67188 21.54649 31.79726 # P GBGR
2 26.18750 21.15291 31.22209 # P GRGB
3 35.95312 29.93572 41.97053 # P GBRG
4 34.60938 28.76076 40.45799 # P GRBG
5 19.93750 13.79119 26.08381 # C GBGR
6 16.39062 10.73761 22.04364 # C GRGB
7 23.35938 16.63430 30.08445 # C GBRG
8 22.79688 16.17650 29.41725 # C GRBG



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
                           c("1" = "GBGR-P", "2"="GRGB-P", "3" = "GBRG-P", 
                             "4" = "GRBG-P", "5" = "GBGR-C", "6" = "GRGB-C",
                             "7" = "GBRG-C", "8" = "GRBG-C"))
F_tall$q.type.cat = revalue(x = as.factor(F_tall$q.type.cat), 
                            c("1" = "Perceptual Question", "2"="Causal Question"))
F_tall$condition.2 = rep(c(1:4), times = 128)
F_tall$condition.2 = revalue(x = as.factor(F_tall$condition.2), 
                             c("1" = "GBGR", "2"="GRGB", "3" = "GBRG", 
                               "4" = "GRBG"))
F_tall$q.type.cat.2 = rep(c(0,1), each = 4, times = 64)
F_tall$q.type.cat.2 = revalue(x = as.factor(F_tall$q.type.cat.2), 
                              c("0" = "Perceptual Question", "1"="Causal Question"))


# create new measure column, 'measure.2', in which each value is subtracted from 100 to
# obtain "inconsistency ratings"
F_tall$measure.2 = (100-F_tall$measure)


#################################### OLD OMNIBUS FIGURE ####################################################
# OMNIBUS ANALYSIS FIGURE
condition_barplot = ggplot(F_tall, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  #facet_wrap(~q.type.cat, scales="free") + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  # PERCEPTUAL SIGNIFICANCE LINES
  geom_signif(comparisons = list(c("GBGR-P", "GBRG-P")), annotations=c("p < .0001"), y_position = 46.5, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GBGR-P", "GRBG-P")), annotations=c("p < .0001"), y_position = 48, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GRGB-P", "GBRG-P")), annotations=c("p < .0001"), y_position = 45, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GRGB-P", "GRBG-P")), annotations=c("p < .0001"), y_position = 43, tip_length = 0.00375) +
  # CAUSAL SIGNIFICANCE LINES
  geom_signif(comparisons = list(c("GRGB-C", "GBRG-C")), annotations=c("p < .05"), y_position = 34, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GRGB-C", "GRBG-C")), annotations=c("p < .07"), y_position = 32, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GBGR-C", "GRGB-C")), annotations=c("p < .0001"), y_position = 28, tip_length = 0.00375) +
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 50)) +
  theme_classic() +
  scale_fill_manual(values=c("#000000", "#999999")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials") # change the main x-axis label
############################################################################################################


#################################### NEW OMNIBUS FIGURE ####################################################
condition_barplot = ggplot(F_tall, aes(q.type.cat.2, measure.2, fill = condition.2)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  geom_signif(annotations = c("p < .0001","p < .0001","p < .0001","p < .0001"),
              y_position = c(49,47,45,43), xmin=c(.6,.875,.6,.875), 
              xmax=c(1.35,1.35,1.1,1.1), 
              tip_length = 0.00375) +
  geom_signif(annotations = c("p < .07","p < .05"),
              y_position = c(34,32), xmin=c(2.275,2.1), 
              xmax=c(1.875,1.875), 
              tip_length = 0.00375) +
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 51)) +
  theme_classic() +
  scale_fill_manual(values = c("white", "gray81", "gray38", "black")) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")   
############################################################################################################



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
chi.data = matrix(c(18,15,23,10,4,3,12,9,7,25,0,2,34,55),2)
dimnames(chi.data) = list(c("Perceptual", "Causal"), c("Markov", "Independence", "Temporal",
                                                       "Other","All Familiar","All Inconsistent",
                                                       "All Associative"))
## COUNTS
# Percep: M=18, I=23, T=4, O=12, F=7, AI=0, A = 34
# Causal: M=15, I=10, T=3, O=9, F=25, AI=2, A = 55


  # LEGEND: 
            # M = Markov; I = Independence; 
            # T = Temporal Precedence; O = Other; 
            # F = ALL Familiar, AI = All Inconsistent
            # A = Associative (sum of all associative options)

# check structure of data
str(chi.data)


# run chisq.test() on chi table
chi.test = chisq.test(chi.data, simulate.p.value = TRUE) # used 'simulate.p.value' 
                                                         # because some cell counts are small

## run post-hoc tests (where you average over columns) ##
for(i in 1:6) print(chisq.test(chi.data[,i])) # this gives you all column-wise comparisons 
                                              # (i.e., M v. I, M v. O, I v. O)

# the chi data
              Markov Independence Temporal Other All Familiar All Inconsistent All Associative
Perceptual     18           23        4    12            7                0              34
Causal         15           10        3     9           25                2              55



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

                  ##
## CAUSAL CONDITION POST-HOC TESTS: ##
                  ##

# Independence comparisons
# M v I
binom_func(15,10)
# I v T
binom_func(10,3)
# I v O
binom_func(10,9)
# I v F
binom_func(10,25)



# Markov comparisons
# M v T
binom_func(15,3)
# M v O
binom_func(15,9)

# M v IC
binom_func(15,2)

# M v AA
binom_func(15,55)



# Collapsing over question type
binom_func(15,55)





