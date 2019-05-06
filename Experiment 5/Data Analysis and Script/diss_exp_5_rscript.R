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
D$med.split.age = rep(0,length(D$age))
for(i in 1:length(D$age)){
  compute = ifelse(D$age[i]<=19.02,0,1)
  D$med.split.age[i] = compute
}

  
# add median split for num.hab
D$med.split.num.hab = rep(0, length(D$num.hab))
for(i in 1:length(D$num.hab)){
  calc = ifelse(D$num.hab[i]<=5,0,1)
  D$med.split.num.hab[i] = calc
}
  
# remove unneeded columns
D$Partic..= NULL
names(D)

# reorder columns
D = as.data.frame(D[,c(1,2,6,22,3,4,
                       5,7,23,8,9:21)])
names(D)


# convert data from "wide" format to "tall" format
D_tall = reshape(D, varying = c(18:21,23), v.names = "measure", 
                 timevar = "test.trial.level", idvar = "ID", 
                 direction = "long")

D_tall = D_tall[order(D_tall$ID),]



# FACTORIZE 'SEX'
D_tall$sex = revalue(x = as.factor(D_tall$sex), 
                     c("0" = "F", "1"="M"))

# FACTORIZE 'TEST.STIM.ORDER'
D_tall$test.stim.order = revalue(x = as.factor(D_tall$test.stim.order),
                                 c("0" = "1234", "1" = "2413"))


# FACTORIZE 'MED.SPLIT.NUM.HAB'
D_tall$med.split.num.hab = revalue(x = as.factor(D_tall$med.split.num.hab),
                                 c("0" = "Fast Habituators", "1" = "Slow Habituators"))


# FACTORIZE 'GROUP'
D_tall$group = revalue(x = as.factor(D_tall$group),
                                   c("0" = "Red", "1" = "Blue"))

# FACTORIZE 'HAB.STIM.ORDER'
D_tall$hab.stim.order = revalue(x = as.factor(D_tall$hab.stim.order), 
                                 c("0" = "12", "1"="21", "2"="34", "3"="43"))


# FACTORIZE 'TEST.TRIAL.LEVEL'
D_tall$test.trial.level = revalue(x = as.factor(D_tall$test.trial.level), 
                                c("1" = "gbgr", "2"="grgb", 
                                  "3"="gbrg", "4"="grbg",
                                  "5"="Posttest"))


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
shapiro.ps = rep(0,2)
for(i in c("Red","Blue")) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$group==i])
  shapiro.ps[i] = shap.calc$p.value
}


## HOMOSKEDASTICITY CHECKS
# plot the boxplots

# homoskedasticity check for "red" group
boxplot(D_tall$measure[D_tall$test.stim.order==0]~D_tall$group[D_tall$group==c(1:4) & D_tall$test.stim.order==0])

#homoskedasticity check for the "red" and "blue" groups
boxplot(D_tall$measure[D_tall$group=="Red"]~D_tall$test.trial.level[D_tall$group=="Red"])
boxplot(D_tall$measure[D_tall$group=="Blue"]~D_tall$test.trial.level[D_tall$group=="Blue"])

# formal test of equal variance for both the 'red' and 'blue' groups
# red group
leveneTest(D_tall$measure[D_tall$group=="Red"]~as.factor(D_tall$test.trial.level[D_tall$group=="Red"]))

# blue group
leveneTest(D_tall$measure[D_tall$group=="Blue"]~as.factor(D_tall$test.trial.level[D_tall$group=="Blue"]))



## ASSUMPTION CHECK NOTES ##
# Given that there is evidence of non-normality, but not heteroskedasticity,
# non-parametric boostrapping and permuation
# testing will be used to estimate confidence intervals and for hypothesis checking



##############################
#### PRELIMINARY ANALYSIS ####
##############################
# analysis to determine whether there is an effect of sex, test trial order, or group on 
# looking time (in s)
lme.fit.prelim = lme(measure~(test.trial.level+group+test.stim.order)^3,random=~1|ID, 
                     na.action=na.exclude, data=D_tall)

lme.fit.prelim = lme(measure~sex+test.stim.order+group+sex:test.stim.order:group, 
                     random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)


########################
# HABITUATION ANALYSIS #
########################
# average and range age
mean(D_tall$age)
range(D_tall$age)

# average number of trials to habituate
mean(D_tall$num.hab)

# num of males
sum(D$sex)


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
dif.num.hab.mean = 6.28125  + 1.96*c(-0.1883883, 0.1883883)

# bootstrapped to obtain CI for mean(D_tall$ttl.hab)
set.seed(2018)
ttl.hab.mean = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d$ttl.hab) 
  return(dif.1)
}

dif.ttl.hab = boot(D_tall, ttl.hab.mean, R=5000) 
dif.ttl.hab
dif.ttl.hab.mean = 84.22813  + 1.96*c(-3.848969, 3.848969)



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
P_tall = reshape(D, varying = c(22:23), v.names = "measure", 
                   timevar = "test.trial.level", idvar = "ID", 
                   direction = "long")


P_tall = P_tall[order(P_tall$ID),]


# set appropriate factor variables in "tall" data
P_tall$sex = as.factor(P_tall$sex)
P_tall$group = as.factor(P_tall$group)
P_tall$age = as.factor(P_tall$age)
P_tall$med.split.age = as.factor(P_tall$med.split.age)
P_tall$hab.stim.order = as.factor(P_tall$hab.stim.order)
P_tall$test.stim.order = as.factor(P_tall$test.stim.order)
P_tall$test.trial.level = as.factor(P_tall$test.trial.level)

P_tall$hab.min.2 = NULL
P_tall$hab.min.3 = NULL


lme.fit.prepost = lme(measure~(test.trial.level+group+sex+test.stim.order+test.trial.level)^4,random=~1|ID, 
                     na.action=na.exclude, data=P_tall)
anova.lme(lme.fit.prepost)

# pretest mean and SD
mean(P_tall$measure[P_tall$test.trial.level=="1"])
sd(P_tall$measure[P_tall$test.trial.level=="1"])

# posttest mean and SD
mean(P_tall$measure[P_tall$test.trial.level=="2"], na.rm = TRUE)
sd(P_tall$measure[P_tall$test.trial.level=="2"], na.rm = TRUE)


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
########################
# MAIN ANCOVA ANALYSIS #
########################
########################

## age as a covariate ##
lme.fit.main = lme(measure~test.trial.level, 
                     random=~1|ID, na.action=na.exclude,
                   data=D_tall)
anova.lme(lme.fit.main)


# for() loop to get means of trials
means_vec = rep(0,5)
for(i in c("gbgr","grgb","gbrg","grbg",
           "Posttest")){
  calc = mean(D_tall$measure[D_tall$test.trial.level==i], 
              na.rm=TRUE)
  means_vec[i] = calc
}
means_vec


# means: 10.865625  7.115625  10.078125  6.775000 19.333548


# nest for() loop to get differences between means
D_tall_2 = reshape(D, varying = c(18:21,23), v.names = "measure", 
                 timevar = "test.trial.level", idvar = "ID", 
                 direction = "long")


D_tall_2 = D_tall_2[order(D_tall_2$ID),]

k=1
comb_vec = rep(0,5)
for(i in 1:5){
  for(j in 1:5){
    calc = (mean(D_tall_2$measure[D_tall_2$test.trial.level==i],na.rm = TRUE)-
              mean(D_tall_2$measure[D_tall_2$test.trial.level==j],na.rm = TRUE))
    comb_vec[k] = calc
    k = k+1
  }
}
comb_vec

# DIFFERENCES BETWEEN PAIRS OF MEANS #
[1]   0.000000   3.750000   0.787500   4.090625  -8.467923  -3.750000   0.000000  -2.962500   0.340625 -12.217923  -0.787500   2.962500   0.000000   3.303125  -9.255423  -4.090625  -0.340625  -3.303125
[19]   0.000000 -12.558548   8.467923  12.217923   9.255423  12.558548   0.000000

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
    dif.1 =  mean(x, data=D_tall_2, na.rm=TRUE) 
    return(dif.1)
  }
  
  boot_main = boot(D_tall_2, boot_func, R=5000) 
  boot_mean[i,] = c(boot_main$t0, boot_main$t0  + 1.96*-sd(boot_main$t), 
                    boot_main$t0  + 1.96*sd(boot_main$t))
}
boot_mean

######################
  # BOOT TEST DF #
######################
V1       V2        V3
1 10.865625  7.498880 14.23237
2  7.115625  4.872700  9.35855
3 10.078125  6.806933 13.34932
4  6.775000  3.976870  9.57313
5 19.333548 15.567131 23.09997




##############################################
## FOLLOW UP PLANNED COMPARISONS PERM TESTS ##
##############################################
# GENERAL PERMUTATION FUNCTION
main_perm_func = function(a,z){
  b = rep(0,5000) 
  for(i in 1:length(b)){
    x = factor(D_tall$test.trial.level, levels=c(a,z)) 
    y = sample(D_tall$measure, replace=TRUE) 
    lm_1 = lm(y ~ x, data=D_tall)
    b[i] = coef(lm_1)[2]
  }
  bb_dif = mean(D_tall$measure[D_tall$test.trial.level==a], na.rm=TRUE)-
    mean(D_tall$measure[D_tall$test.trial.level==z], na.rm=TRUE)
  c((sum(abs(b) > bb_dif)/length(b)),(sum(b > bb_dif)/length(b)),
    sum((abs(b) < bb_dif)/length(b)),(sum(b < bb_dif)/length(b)),
    mean(D_tall$measure[D_tall$test.trial.level==a], na.rm=TRUE),
    mean(D_tall$measure[D_tall$test.trial.level==z], na.rm=TRUE),
    bb_dif)
}

################################
# PERM TESTS: MAIN TEST EVENTS #
################################
# levels of D_tall$test.trial
levels(D_tall$test.trial)
[1] "gbgr"     "grgb"     "gbrg"     "grbg"     "Posttest"

# 
# GBGR V GRGB
main_perm_func("gbgr","grgb")

# GBGR V GBRG
main_perm_func("gbgr","gbrg")

# GBGR V GRBG
main_perm_func("gbgr","grbg")

# GRGB V GBRG
main_perm_func("grgb","gbrg")

# GRGB V GRBG
main_perm_func("grgb","grbg")

# GBRG V GRBG
main_perm_func("gbrg","grbg")


###########################################################
# PERM TESTS: POST TEST EVENT COMPARED TO ALL TEST EVENTS #
###########################################################
# POSTTEST V GBGR
main_perm_func("Posttest","gbgr")

# POSTTEST V GRGB
main_perm_func("Posttest","grgb")

# POSTTEST V GBRG
main_perm_func("Posttest","gbrg")

# POSTTEST V GRBG
main_perm_func("Posttest","grbg")


##################################################################
# BAYES FACTOR (FUNCTION) FOR MARGIANLLY SIGNIFICANT DIFFERENCES #
##################################################################
bayes_factor_func = function(a,b,x){
  modified_df = subset(D_tall, ! test.trial.level %in% c(a,b,x))
  lm.null = lme(measure~1, random=~1|ID, data=modified_df)
  lm.alt = lme(measure~test.trial.level, random=~1|ID, data=modified_df)
  
  null.bic = BIC(lm.null)
  alt.bic = BIC(lm.alt)
  
  # compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)
  
  BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
  BF10 = 1/BF01
  BF10
}

# gbgr vs grgb
bayes_factor_func(a="gbrg",b="grbg",x="Posttest")
# 3.531591

# grgb vs grbg
bayes_factor_func(a="grgb",b="gbrg",x="Posttest")
# 3.739883

# gbrg vs grbg
bayes_factor_func(a="gbgr",b="grgb",x="Posttest")
# 2.487306

# gbgr vs grgb
bayes_factor_func(a="grgb",b="grbg",x="Posttest")
#0.6499908

# grbg vs grgb
bayes_factor_func(a="gbgr",b="gbrg",x="Posttest")
# 0.6104696

# gbgr vs grbg
bayes_factor_func(a="grgb",b="gbrg",x="Posttest")
# 0.6104696

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
D_tall$med.split.age = revalue(x = as.factor(D_tall$med.split.age), 
                                  c("0" = "Younger (<=19.01)", "1" = "Older (>19.01)"))

# OMNIBUS ANALYSIS FIGURE
condition_barplot = ggplot(D_tall, aes(test.trial.level, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
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
main.age.fit = lme(measure~(test.trial.level+med.split.age)^2, 
                   random=~1|ID, na.action=na.exclude,
                   data=D_tall)

anova.lme(main.age.fit)

# create separate DFs for both the younger and older infants
U = subset(D_tall, ! test.trial.level %in% c("Posttest"))
U = U[order(U$med.split.age),]


younger_df = U[c(1:68),]
older_df = U[c(69:128),]  


# Omnibus Younger and Older Plot
condition_barplot = ggplot(U, aes(med.split.age, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
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
  theme(legend.box.background = element_rect(), legend.box.margin = margin(3, 3, 3, 3)) + # this creates a square around legends
  theme(legend.text = element_text(size = 12)) + # this sets the size of the legend text 
  theme(legend.title=element_blank()) + # this removes the automatic legend title
  labs(x = "Test trials")



# Younger Infants
condition_barplot = ggplot(younger_df, aes(test.trial.level, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("white", "gray81", "gray38", "gray20", "black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  labs(x = "Test trials")
  

# Older Infants
condition_barplot = ggplot(older_df, aes(test.trial.level, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("white", "gray81", "gray38", "gray20", "black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
labs(x = "Test trials")




###################################################################
## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR OLDER INFANTS ##
###################################################################
# install and load post-hoc chi-square test package (called follow up binomial tests)

# create contigency table with counts
chi.younger.data = matrix(c(6,3,8),1)
dimnames(chi.younger.data) = list(c("Younger Infant"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Younger: M=6, P=3,  O = 8


# the chi data
                 Markov Perceptual Other
Younger Infant      6          3     8

# run chisq.test() on chi table
chi.younger.test = chisq.test(chi.younger.data, simulate.p.value = TRUE) 
chi.younger.test

# YOUNGER INFANTS INDIVIDUAL DIFFERENCES SUMMARY
# Younger infants did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.



###################################################################
## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR OLDER INFANTS ##
###################################################################
# install and load post-hoc chi-square test package (called follow up binomial tests)

# create contigency table with counts
chi.older.data = matrix(c(3,4,8),1)
dimnames(chi.older.data) = list(c("Younger Infant"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Older: M=3, P=4,  O = 8


# the chi data
                Markov Perceptual Other
Older Infant      3          4     8

# run chisq.test() on chi table
chi.older.test = chisq.test(chi.older.data, simulate.p.value = TRUE) 
chi.older.test

# OLDER INFANTS INDIVIDUAL DIFFERENCES SUMMARY
# Older infants, like younger infants, did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.





##########################################################
## INDIVIDUAL DIFFERENCES FOR FAST AND SLOW HABITUATORS ##
##########################################################
main.hab.rate.fit = lme(measure~(test.trial.level+med.split.num.hab)^2, 
                   random=~1|ID, na.action=na.exclude,
                   data=D_tall)

anova.lme(main.hab.rate.fit)



# create separate DFs for both the younger and older infants
G = subset(D_tall, ! test.trial.level %in% c("Posttest"))
G = G[order(G$med.split.num.hab),]


fast_hab_df = G[c(1:60),]
slow_hab_df = G[c(61:128),] 

# HABITUATION RATE PERMUTATION FUNCTION
hab_rate_perm_func = function(a,z,n){
  b = rep(0,5000) 
  for(i in 1:length(b)){
    x = factor(G$test.trial.level, levels=c(a,z))
    y = sample(G$measure, replace=TRUE) 
    lm_1 = lm(y[G$med.split.num.hab==n] ~ x[G$med.split.num.hab==n], data=G)
    b[i] = coef(lm_1)[2]
  }
  bb_dif = mean(G$measure[G$test.trial.level==a & G$med.split.num.hab==n], na.rm=TRUE)-
    mean(G$measure[G$test.trial.level==z & G$med.split.num.hab==n], na.rm=TRUE)
  c((sum(abs(b) > bb_dif)/length(b)),(sum(b > bb_dif)/length(b)),
    sum((abs(b) < bb_dif)/length(b)),(sum(b < bb_dif)/length(b)),
    mean(G$measure[G$test.trial.level==a & G$med.split.num.hab==n], na.rm=TRUE),
    mean(G$measure[G$test.trial.level==z & G$med.split.num.hab==n], na.rm=TRUE),
    bb_dif)
}


# Fast Habituators Planned Comparisons
# gbgr v grgb
hab_rate_perm_func("gbgr","grgb", "Fast Habituators")
# gbgr v gbrg
hab_rate_perm_func("gbgr","gbrg", "Fast Habituators")
# gbgr v grbg
hab_rate_perm_func("gbgr","grbg", "Fast Habituators")
# grgb v gbrg
hab_rate_perm_func("grgb","gbrg", "Fast Habituators")
# grgb v grbg
hab_rate_perm_func("grgb","grbg", "Fast Habituators")
# gbrg v grbg
hab_rate_perm_func("gbrg","grbg", "Fast Habituators")


# Slow Habituators Planned Comparisons
# gbgr v grgb
hab_rate_perm_func("gbgr","grgb", "Slow Habituators")
# gbgr v gbrg
hab_rate_perm_func("gbgr","gbrg", "Slow Habituators")
# gbgr v grbg
hab_rate_perm_func("gbgr","grbg", "Slow Habituators")
# grgb v gbrg
hab_rate_perm_func("grgb","gbrg", "Slow Habituators")
# grgb v grbg
hab_rate_perm_func("grgb","grbg", "Slow Habituators")
# gbrg v grbg
hab_rate_perm_func("gbrg","grbg", "Slow Habituators")



# Fast and Slow Habituators Bootstrapped Means
hab_boot_mean_func = function(x){
  hab_boot_mean = as.data.frame(matrix(NA, nrow=4, ncol=3, byrow=TRUE))
  for(i in 1:nrow(hab_boot_mean)){ # want number of iterations to equal number of rows, especially because we're filling in by row
    set.seed(2018)
    boot_func = function(data,b,formula, p){ 
      d= data[b,] 
      x = d$measure[d$test.trial.level==i & d$med.split.num.hab==x]
      dif.1 =  mean(x, data=D_tall_2, na.rm=TRUE) 
      return(dif.1)
    }
    
    boot_main = boot(D_tall_2, boot_func, R=5000) 
    hab_boot_mean[i,] = c(boot_main$t0, boot_main$t0  + 1.96*-sd(boot_main$t), 
                          boot_main$t0  + 1.96*sd(boot_main$t))
  }
  hab_boot_mean
}


# SLOW HABITUATORS BOOTSTRAPPED CIs
hab_boot_mean_func(1)
V1       V2        V3
1 13.794118 8.425743 19.162492
2  5.964706 2.972546  8.956866
3 10.811765 6.052085 15.571444
4  4.523529 2.593708  6.453350

# FAST HABITUATORS BOOTSTRAPPED CIs
hab_boot_mean_func(0)
V1       V2       V3
1 7.546667 4.339970 10.75336
2 8.420000 5.032977 11.80702
3 9.246667 4.741503 13.75183
4 9.326667 3.841989 14.81134

# OMNIBUS HABITUATORS GRAPH
condition_barplot = ggplot(G, aes(med.split.num.hab, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  #facet_wrap(~q.type.cat, scales="free") + # create as many separate graphs as there are conditions 
  ylab("Looking time (s)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 25)) +
  theme_classic() +
  scale_fill_manual(values = c("white", "gray81", "gray38", "gray20", "black")) + # changes colors of individual bars
  theme(legend.box.background = element_rect(), legend.box.margin = margin(3, 3, 3, 3)) + # this creates a square around legends
  theme(legend.text = element_text(size = 12)) + # this sets the size of the legend text 
  theme(legend.title=element_blank()) + # this removes the automatic legend title
  labs(x = "Test trials") 

# SLOW HABITUATORS
condition_barplot = ggplot(slow_hab_df, aes(test.trial.level, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("white", "gray81", "gray38", "gray20", "black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
labs(x = "Test trials")



# FAST HABITUATORS
condition_barplot = ggplot(fast_hab_df, aes(test.trial.level, measure, fill = test.trial.level)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  facet_wrap(~ID) + # create as many separate graphs as there are conditions 
  ylab("Looking Time (s)") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 80)) +
  theme(strip.background =element_rect(fill='black')) +
  theme(strip.text = element_text(colour = 'white', size = 12, face = "bold")) +
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("white", "gray81", "gray38", "gray20", "black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
labs(x = "Test trials")


####################################################################
## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR SLOW HABITUATORS ##
####################################################################
# install and load post-hoc chi-square test package (called follow up binomial tests)

# create contigency table with counts
chi.slow.data = matrix(c(4,2,11),1)
dimnames(chi.slow.data) = list(c("Slow Habituators"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Older: M=4, P=2,  O = 11


# the chi data
                    Markov Perceptual Other
Slow Habituators      4          2     11

# run chisq.test() on chi table
chi.slow.test = chisq.test(chi.slow.data, simulate.p.value = TRUE)


# Follow-up binomial tests
binom.test(4,15,p=1/2) # Markov v Other
binom.test(2,13,p=1/2) # Perceptual v other
binom.test(4,6,p=1/2) # Markov v Perceptual

# SLOW HABITUATORS INDIVIDUAL DIFFERENCES SUMMARY
# Older infants, like younger infants, did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.



####################################################################
## CHI SQUARE TEST ON INDIVIDUAL DIFFERENCES FOR FAST HABITUATORS ##
####################################################################
# install and load post-hoc chi-square test package (called follow up binomial tests)

# create contigency table with counts
chi.fast.data = matrix(c(3,3,9),1)
dimnames(chi.fast.data) = list(c("Fast Habituators"), c("Markov", "Perceptual", "Other"))
## COUNTS
# Older: M=5, P=2,  O = 6


# the chi data
Markov Perceptual Other
Fast Habituators      3          3     9

# run chisq.test() on chi table
chi.fast.test = chisq.test(chi.fast.data, simulate.p.value = TRUE) 
chi.fast.test


# Follow-up binomial tests
binom.test(3,12,p=1/2) # Markov v Other
binom.test(3,12,p=1/2) # Perceptual v other
binom.test(3,6,p=1/2) # Markov v Perceptual

# FAST HABITUATORS INDIVIDUAL DIFFERENCES SUMMARY
# Older infants, like younger infants, did not differ in the extent to which they 
# were classified as Markov processors, Perceptual processors, 
# Other processors.