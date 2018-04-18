########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 3 SCRIPT     #############
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
D$X = NULL
D = D[c(1:32),]

# reorder columns
D = as.data.frame(D[,c(1,2,3,5,7,9,4,6,8,10,11)])


# reshape the data
D_tall = reshape(D, varying = 3:10, v.names = "measure", 
                 timevar = "condition", idvar = "ID", 
                 new.row.names = 1:256, direction = "long")

# order data
D_tall = D_tall[order(D_tall$ID),]

# add q.type.cat column
D_tall$q.type.cat = as.factor(rep(c(1:2), each = 4, times = 32))

# names D_tall
"ID"         "exp"        "q.type"     "condition"  "measure"    "q.type.cat"


# set appropriate factor variables in "tall" data
D_tall$condition = as.factor(D_tall$condition)
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
# Given that there is evidence of heteroskedasticity and non-normality,
# non-parametric boostrapping and permuation
# testing will be used to estimate confidence intervals and for hypothesis checking



##############################
#### PRELIMINARY ANALYSIS ####
##############################
# analysis to determine effect of question type and 
lme.fit.prelim = lme(measure~q.type.cat+q.type+q.type.cat:q.type, random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)




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
    x = d[d$condition==i,7]
    dif.1 =  mean(x, data=D_tall) 
    return(dif.1)
  }
  
  GBGR_P_boot = boot(D_tall, boot_func, R=4000) 
  boot_mean[i,] = c(GBGR_P_boot$t0, GBGR_P_boot$t0  + 1.96*-sd(GBGR_P_boot$t), 
                    GBGR_P_boot$t0  + 1.96*sd(GBGR_P_boot$t))
}


#####################################
## BOOT TESTS: PERCEPTUAL & CAUSAL ##
#####################################
V1         V2        V3
1  2.81250  0.4474658  5.177534 # GFCF-P
2 20.46875 14.5977934 26.339707 # GFCN-P
3 38.81250 30.6207063 47.004294 # GNCF-P
4 33.93750 26.5511353 41.323865 # GNCN-P
5  3.81250  0.5068676  7.118132 # GFCF-C
6 19.06250 12.4821339 25.642866 # GFCN-C
7 33.65625 26.3920905 40.920410 # GNCF-C
8 29.62500 23.3216403 35.928360 # GNCN-C



############################
## PERM TESTS: PERCEPTUAL ##
############################
# GFCF-P v. GFCN-P
perm_func(1,2)

# GFCF-P v. GNCF-P
perm_func(1,3)

# GFCF-P v. GNCN-P
perm_func(1,4)

# GFCN-P v. GNCF-P
perm_func(2,3)

# GFCN-P V. GNCN-P
perm_func(2,4)

# GNCF-P V. GNCN-P
perm_func(3,4)


############################
## PERM TESTS: CAUSAL ##
############################
# GFCF-C v. GFCN-C
perm_func(5,6)

# GFCF-C v. GNCF-C
perm_func(5,7)

# GFCF-C v. GNCN-C
perm_func(5,8)

# GFCN-C v. GNCF-C
perm_func(6,7)

# GFCN-C V. GNCN-C
perm_func(6,8)

# GNCF-C V. GNCN-C
perm_func(7,8)


#####################################################
# Bayes Factor to compare GNCF v GNCN CAUSAL EVENTS #
#####################################################
sub.bayes = subset(D_tall, ! condition %in% c(1:5,7))
# define the null and alternative models #
lm.null = lme(measure~1, random=~1|ID, data=sub.bayes)
lm.alt = lm(measure~as.factor(condition), data=sub.bayes)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01


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
                           c("1" = "GFCF-P", "2"="GFCN-P", "3" = "GNCF-P", 
                             "4" = "GNCN-P", "5" = "GFCF-C", "6" = "GFCN-C",
                             "7" = "GNCF-C", "8" = "GNCN-C"))
F_tall$q.type.cat = revalue(x = as.factor(F_tall$q.type.cat), 
                            c("1" = "Perceptual Question", "2"="Causal Question"))


# OMNIBUS ANALYSIS FIGURE
condition_barplot = ggplot(F_tall, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  #facet_wrap(~q.type.cat, scales="free") + # create as many separate graphs as there are conditions 
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  # PERCEPTUAL SIGNIFICANCE LINES
  geom_signif(comparisons = list(c("GFCF-P", "GFCN-P")), annotations=c("p < .001"), y_position = 30, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCF-P", "GNCF-P")), annotations=c("p < .00001"), y_position = 57, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCF-P", "GNCN-P")), annotations=c("p < .00001"), y_position = 54, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCN-P", "GNCF-P")), annotations=c("p < .005"), y_position = 51, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCN-P", "GNCN-P")), annotations=c("p < .0275"), y_position = 48, tip_length = 0.00375) +
  # CAUSAL SIGNIFICANCE LINES
  geom_signif(comparisons = list(c("GFCF-C", "GFCN-C")), annotations=c("p < .005"), y_position = 30, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCF-C", "GNCF-C")), annotations=c("p < .00001"), y_position = 52, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCF-C", "GNCN-C")), annotations=c("p < .00001"), y_position = 49, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCN-C", "GNCF-C")), annotations=c("p < .01"), y_position = 45, tip_length = 0.00375) +
  geom_signif(comparisons = list(c("GFCN-C", "GNCN-C")), annotations=c("p < .05"), y_position = 42, tip_length = 0.00375) +
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 75)) +
  theme_classic() +
  scale_fill_manual(values=c("#000000", "#999999")) +
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
condition_barplot = ggplot(F_tall, aes(condition, measure.2, fill = q.type.cat)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
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
chi.data = matrix(c(2,2,19,20,1,5,2,1,1,0,7,4),2)
dimnames(chi.data) = list(c("Perceptual", "Causal"), c("Markov", "Perceptual", "Causal",
                                                       "Last_3_Novel","All Familiar","Other"))
## COUNTS
# Percep: M=2, P=19, C=1, l3N=2, AF = 1, O = 7
# Causal: M=2, P=20, C=5, l3N=1, AF = 0, O = 4

# the chi data
              Markov Perceptual Causal Last_3_Novel All Familiar Other
Perceptual      2         19      1            2            1     7
Causal          2         20      5            1            0     4


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
# P v M
binom_func(2,21)
0.00006604195

# P v C
binom_func(1,21)
0.00001096725

# P v L3N
binom_func(2,21)
0.00006604195

# P v AF
binom_func(1,21)
0.00001096725

# P v O
binom_func(7,21)
0.0275


# M v L3N
binom_func(2,4)
0.6875

# M v AF 
binom_func(1,3)
0.625

# M v 0
binom_func(2,9)
0.06542969





                      ##
## CAUSAL CONDITION POST-HOC TESTS: ##
                      ##
# P v M
binom_func(2,22)
0.000035882

# P v C
binom_func(5,25)
0.0003249142

# P v L3N
binom_func(1,21)
0.00001096725

# P v AF
binom_func(0,21)
0.0000009536743

# P v O
binom_func(4,24)
0.0001799911

# M v L3N
binom_func(1,3)
0.625

# M v AF 
binom_func(0,2)
0.5

# M v 0
binom_func(2,6)
0.2890625
