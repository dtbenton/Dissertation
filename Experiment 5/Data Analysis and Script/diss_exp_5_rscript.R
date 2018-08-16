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
library(ez)
options(scipen=9999)

# load data
D = read.csv(file.choose(), header = TRUE)


# add median split age column
D$med.split.age = rep(0,25)
for(i in 1:25){
  compute = ifelse(D$age[i]<=18.06,0,1)
  D$med.split.age[i] = compute
}

  
# add median split for num.hab
D$med.split.num.hab = rep(0, 25)
for(i in 1:25){
  calc = ifelse(D$num.hab[i]<=5,0,1)
  D$med.split.num.hab[i] = calc
}
  
# reorder columns
D = as.data.frame(D[,c(1:6,22,7:23)])

# delete redundant med.split.age column
D$med.split.age.1 = NULL

# convert data from "wide" format to "tall" format
D_tall = reshape(D, varying = c(17:20,22), v.names = "measure", 
                 timevar = "test.trial.level", idvar = "ID", 
                 direction = "long")


D_tall = D_tall[order(D_tall$ID),]

# set appropriate factor variables in "tall" data
D_tall$sex = as.factor(D_tall$sex)
D_tall$group = as.factor(D_tall$group)
D_tall$med.split.age = as.factor(D_tall$med.split.age)
D_tall$hab.stim.order = as.factor(D_tall$hab.stim.order)
D_tall$test.stim.order = as.factor(D_tall$test.stim.order)
D_tall$test.trial.level = as.factor(D_tall$test.trial.level)

# reorder columns of D_tall
D_tall = as.data.frame(D_tall[,c(1:2,4:18,3,19)])



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


## HOMOSKEDASTICITY CHECKS
# plot the boxplots

# homoskedasticity check for "red" group
boxplot(D_tall$measure[D_tall$test.stim.order==0]~D_tall$group[D_tall$group==c(1:4) & D_tall$test.stim.order==0])

#homoskedasticity check for "blue" group
boxplot(D_tall$measure[D_tall$test.stim.order==1]~D_tall$group[D_tall$group==c(1:4) & D_tall$test.stim.order==1])


# formal test of equal variance for both the 'red' and 'blue' groups
# red group
leveneTest(D_tall$measure[D_tall$group==0]~as.factor(D_tall$test.trial.level[D_tall$group==0]))

# blue group
leveneTest(D_tall$measure[D_tall$group==1]~as.factor(D_tall$test.trial.level[D_tall$group==1]))



## ASSUMPTION CHECK NOTES ##
# Given that there is evidence of non-normality, but not heteroskedasticity,
# non-parametric boostrapping and permuation
# testing will be used to estimate confidence intervals and for hypothesis checking



##############################
#### PRELIMINARY ANALYSIS ####
##############################
# analysis to determine whether there is an effect of sex
lme.fit.prelim = lme(measure~sex+test.stim.order+sex:test.stim.order, 
                     random=~1|ID, data=D_tall)

lme.fit.prelim = lme(measure~sex+test.stim.order+group+sex:test.stim.order:group, 
                     random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)


########################
# HABITUATION ANALYSIS #
########################
# average number of trials to habituate
mean(D_tall$num.hab)


# average seconds to habituate
mean(D_tall$ttl.hab)


# bootstrapped to obtain CI for mean(D_tall$num.hab)
set.seed(2018)
num.hab.mean = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d$num.hab) 
  return(dif.1)
}

dif.num.hab = boot(D_tall, num.hab.mean, R=5000) 
dif.num.hab
dif.num.hab.mean = 6.24  + 1.96*c(-0.2196325, 0.2196325)

# bootstrapped to obtain CI for mean(D_tall$ttl.hab)
set.seed(2018)
ttl.hab.mean = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d$ttl.hab) 
  return(dif.1)
}

dif.ttl.hab = boot(D_tall, ttl.hab.mean, R=5000) 
dif.ttl.hab
dif.ttl.hab.mean = 81.352  + 1.96*c(-4.471716, 4.471716)



# EFFECTS OF SEX AND TEST TRIAL ORDER ON TOTAL LOOKING
# TO THE HABITUATION EVENTS
ez.ttl.hab = ezANOVA(D_tall, dv = ttl.hab, between = .(sex, group, test.stim.order), 
                     wid = ID)
print(ez.ttl.hab) 



# EFFECTS OF SEX AND TEST TRIAL ORDER ON NUM OF TRIALS 
# TO REACH HABITUATION CRITERION
ez.num.hab = ezANOVA(D_tall, dv = num.hab, between = .(sex, group, test.stim.order), 
                     wid = ID)
print(ez.num.hab) 

test.stim.order.hab.2 = lme(num.hab~test.stim.order, 
                          random=~1|ID, data=D_tall)
anova.lme(test.stim.order.hab.2)

group.hab.2 = lme(ttl.hab~group, 
                random=~1|ID, data=D_tall)
anova.lme(group.hab.2)


# ANCOVA TO CONTROL FOR AGE WHEN DETERMINING WHETHER
# THERE ARE EFFECTS OF SEX AND TEST TRIAL ORDER ON TOTAL LOOKING
# TO THE HABITUATION EVENTS AND NUMBER OF TRIALS TO REACH 
# HABITUATION CRITERION

# ttl.hab
ancova.ttl.hab = ezANOVA(D_tall, dv = ttl.hab, between = .(sex, group, test.stim.order), 
                     wid = ID,
                     between_covariates=age)
print(ancova.ttl.hab) 


# num.hab
ancova.num.hab = ezANOVA(D_tall, dv = num.hab, between = .(sex, group, test.stim.order),
                         wid = ID,
                         between_covariates=age)
print(ancova.num.hab)




########################
# PRETEST VS. POSTTEST #
########################
P_tall = reshape(D, varying = c(21:22), v.names = "measure", 
                   timevar = "test.trial.level", idvar = "ID", 
                   direction = "long")


P_tall = P_tall[order(P_tall$ID),]


# set appropriate factor variables in "tall" data
P_tall$sex = as.factor(P_tall$sex)
P_tall$group = as.factor(P_tall$group)
P_tall$med.split.age = as.factor(P_tall$med.split.age)
P_tall$hab.stim.order = as.factor(P_tall$hab.stim.order)
P_tall$test.stim.order = as.factor(P_tall$test.stim.order)
P_tall$test.trial.level = as.factor(P_tall$test.trial.level)


prepost = ezANOVA(P_tall, dv = measure, between = .(sex, group, test.stim.order, 
                                                    test.trial.level), 
                         wid = ID,
                         between_covariates=age)
print(prepost) 


## Follow up BF analysis to examine marginally reliable main effect of 'group type' ##
# define the null and alternative models #
lm.null = lme(measure~1, random=~1|ID, data=P_tall)
lm.alt = lm(measure~as.factor(group), data=P_tall)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01


############################
# ANCOVA ASSUMPTION CHECKS #
############################

# check to determine whether the covariate, age, depends on
# the levels of the factor, 'test.trial.level'
age.assump.aov = aov(age~test.trial.level, data=D_tall)
summary(age.assump.aov)

# check homoskedasticity of covariate
ancova.homo.check = lme(measure~test.trial.level+age+test.trial.level:age, 
                        random=~1|ID, data=D_tall)
aov.ancova.homo.check = anova.lme(ancova.homo.check)

# SUMMARY OF ANCOVA ASSUMPTION CHECKS #
# homoskedasticity of covariates and independence of covariate from
# the levels of the factor were not violated. 



########################
# MAIN ANCOVA ANALYSIS #
########################
## age as a covariate ##
ancova.main.fit = ezANOVA(D_tall, dv = measure, within=test.trial.level,
                         wid = ID,
                         between_covariates=age)
print(ancova.main.fit)


## age using a median split ##
ancova.med.split = ezANOVA(D_tall, dv = measure, within = test.trial.level,
                          between = med.split.age,
                          wid = ID)
print(ancova.med.split)


# for() loop to get means of trials
means_vec = rep(0,5)
for(i in 1:length(means_vec)){
  calc= mean(D_tall$measure[D_tall$test.trial.level==i])
  means_vec[i] = calc
}
means_vec

# means: 9.532  7.424  8.956  7.216 20.284


# nest for() loop to get differences between means
k=1
comb_vec = rep(0,25)
for(i in 1:5){
  for(j in 1:5){
    calc = (mean(D_tall$measure[D_tall$test.trial.level==i])-
              mean(D_tall$measure[D_tall$test.trial.level==j]))
    comb_vec[k] = calc
    k = k+1
  }
}


# DIFFERENCES BETWEEN PAIRS OF MEANS #
[1]   0.000   2.108   0.576   2.316 -10.752  -2.108   0.000  -1.532   0.208 -12.860
[11]  -0.576   1.532   0.000   1.740 -11.328  -2.316  -0.208  -1.740   0.000 -13.068
[21]  10.752  12.860  11.328  13.068   0.000

########################
# Global Boot Function #
########################
# bootstrap global function
boot_mean = as.data.frame(matrix(NA, nrow=5, ncol=3, byrow=TRUE))
for(i in 1:nrow(boot_mean)){ # want number of iterations to equal number of rows, especially because we're filling in by row
  set.seed(2018)
  boot_func = function(data,b,formula, p){ 
    d= data[b,] 
    x = d$measure[d$test.trial.level==i]
    dif.1 =  mean(x, data=D_tall) 
    return(dif.1)
  }
  
  boot_main = boot(D_tall, boot_func, R=5000) 
  boot_mean[i,] = c(boot_main$t0, boot_main$t0  + 1.96*-sd(boot_main$t), 
                    boot_main$t0  + 1.96*sd(boot_main$t))
}

######################
  # BOOT TEST DF #
######################
V1        V2       V3
1  9.532  6.516559 12.54744
2  7.424  4.639529 10.20847
3  8.956  5.476158 12.43584
4  7.216  4.088433 10.34357
5 20.284 16.173807 24.39419




##############################################
## FOLLOW UP PLANNED COMPARISONS PERM TESTS ##
##############################################
# GBGR v GRGB
GBGRvGRGB = subset(D_tall, ! test.trial.level %in% c(3:5))
GBGRvGRGB$test.trial.level = factor(GBGRvGRGB$test.trial.level)




# permutation test for B+ mid vs post ratings
b = rep(0,5000) 
for(i in 1:5000){
  y = sample(GBGRvGRGB$measure, replace=TRUE)
  lm_1 = lme(y ~ test.trial.level, random=~1|ID, data=GBGRvGRGB) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~test.trial.level, random=~1|ID, data=GBGRvGRGB)
beta_actual = fixed.effects(lm.fit)[2]

# p value
sum(abs(b) > beta_actual)/5000




########################################################
########################################################
########################################################
#############                              #############
#############            Figures           #############
#############                              #############
########################################################
########################################################
########################################################

## AGE AS COVARIATE ##
# Create 'F_tall' data frame to use for ggplot
F_tall = D_tall

# rename levels of 'condition' and 'q.type.cat' factors
F_tall$test.trial.level = revalue(x = as.factor(F_tall$test.trial.level), 
                           c("1" = "GBGR", "2"="GRGB", "3" = "GBRG", 
                             "4" = "GRBG",
                             "5" = "Posttest"))

F_tall$med.split.age = revalue(x = as.factor(F_tall$med.split.age), 
                                  c("0" = "Younger (<=18.06)", "1" = "Older (>18.06)"))

# OMNIBUS ANALYSIS FIGURE
condition_barplot = ggplot(F_tall, aes(test.trial.level, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  ylab("Looking time (s)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 25)) +
  theme_classic() +
  scale_fill_manual(values = c("white", "gray81", "gray38", "gray20", "black")) + # changes colors of individual bars
  theme(legend.position="none") + # hides legend but still exists in viewport
  labs(x = "Test trials") # change the main x-axis label



##############################################################################
##############################################################################
##############################################################################
#############                                                    #############
#############            INDIVIDUAL DIFFERENCE ANALYSIS          #############
#############                                                    #############
##############################################################################
##############################################################################
##############################################################################



##########################################################
## INDIVIDUAL DIFFERENCES FOR OLDER AND YOUNGER INFANTS ##
##########################################################
ancova.main.age.fit = ezANOVA(D_tall, dv = measure, within=test.trial.level,
                                   between = med.split.age,
                                   wid = ID)
print(ancova.main.age.fit)

# create separate DFs for both the younger and older infants
Z = D_tall
Z = Z[order(Z$med.split.age),]

younger_df = Z[c(1:52),]
older_df = Z[c(53:100),]

# Younger Infants
condition_barplot = ggplot(younger_df, aes(test.trial.level, measure, fill = med.split.age)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
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


## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR YOUNGER INFANTS ##
# install and load post-hoc chi-square test package (called follow up binomial tests)


# create contigency table with counts
chi.younger.data = matrix(c(4,2,7),1)
dimnames(chi.younger.data) = list(c("Younger Infant"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Younger: M=4, P=2,  O = 7


# the chi data
                 Markov Perceptual Other
Younger Infant      4          2     7

# run chisq.test() on chi table
chi.younger.test = chisq.test(chi.younger.data, simulate.p.value = TRUE) 

# YOUNGER INFANTS INDIVIDUAL DIFFERENCES SUMMARY
# Younger infants did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.



# Older Infants
condition_barplot = ggplot(older_df, aes(test.trial.level, measure, fill = med.split.age)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
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

###################################################################
## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR OLDER INFANTS ##
###################################################################
# install and load post-hoc chi-square test package (called follow up binomial tests)

# create contigency table with counts
chi.older.data = matrix(c(4,3,5),1)
dimnames(chi.older.data) = list(c("Younger Infant"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Older: M=4, P=3,  O = 5


# the chi data
                Markov Perceptual Other
Older Infant      4          3     5

# run chisq.test() on chi table
chi.older.test = chisq.test(chi.older.data, simulate.p.value = TRUE) 

# OLDER INFANTS INDIVIDUAL DIFFERENCES SUMMARY
# Older infants, like younger infants, did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.





##########################################################
## INDIVIDUAL DIFFERENCES FOR FAST AND SLOW HABITUATORS ##
##########################################################
ancova.main.hab.rate.fit = ezANOVA(D_tall, dv = measure, within=test.trial.level,
                                   between = med.split.num.hab,
                                   wid = ID)
print(ancova.main.hab.rate.fit)


# create separate DFs for both the younger and older infants
L = D_tall
L = L[order(L$med.split.num.hab),]

# rename levels of 'condition' and 'q.type.cat' factors
L$med.split.num.hab = revalue(x = as.factor(L$med.split.num.hab), 
                              c("0" = "Fast Habituators (<=5)", "1" = "Slow Habituators (>5)"))

fast_habituators_df = L[c(1:52),]
slow_habituators_df = L[c(53:100),]



 
# OMNIBUS HABITUATORS GRAPH
condition_barplot = ggplot(L, aes(med.split.num.hab, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  #facet_wrap(~q.type.cat, scales="free") + # create as many separate graphs as there are conditions 
  ylab("Looking time (s)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 18)) +
  theme_classic() +
  labs(x = "Test trials") 

# SLOW HABITUATORS
condition_barplot = ggplot(slow_habituators_df, aes(test.trial.level, measure, fill = med.split.num.hab)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
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



# FAST HABITUATORS
condition_barplot = ggplot(fast_habituators_df, aes(test.trial.level, measure, fill = med.split.num.hab)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
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


####################################################################
## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR SLOW HABITUATORS ##
####################################################################
# install and load post-hoc chi-square test package (called follow up binomial tests)

# create contigency table with counts
chi.slow.data = matrix(c(5,2,5),1)
dimnames(chi.slow.data) = list(c("Slow Habituators"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Older: M=5, P=2,  O = 5


# the chi data
                    Markov Perceptual Other
Slow Habituators      5          2     6

# run chisq.test() on chi table
chi.slow.test = chisq.test(chi.slow.data, simulate.p.value = TRUE) 

# SLOW HABITUATORS INDIVIDUAL DIFFERENCES SUMMARY
# Older infants, like younger infants, did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.



####################################################################
## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR FAST HABITUATORS ##
####################################################################
# install and load post-hoc chi-square test package (called follow up binomial tests)

# create contigency table with counts
chi.fast.data = matrix(c(5,2,6),1)
dimnames(chi.fast.data) = list(c("Fast Habituators"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Older: M=5, P=2,  O = 6


# the chi data
Markov Perceptual Other
Fast Habituators      5          2     6

# run chisq.test() on chi table
chi.fast.test = chisq.test(chi.fast.data, simulate.p.value = TRUE) 

# FAST HABITUATORS INDIVIDUAL DIFFERENCES SUMMARY
# Older infants, like younger infants, did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.