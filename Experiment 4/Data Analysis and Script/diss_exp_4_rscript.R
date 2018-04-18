########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 4 SCRIPT     #############
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
D = as.data.frame(D[,c(1,2,3,7,5,9,4,8,6,10,11)])



# reshape the data
D_tall = reshape(D, varying = 3:10, v.names = "measure", 
                 timevar = "condition", idvar = "ID", 
                 new.row.names = 1:256, direction = "long")

# order data
D_tall = D_tall[order(D_tall$ID),]

# add q.type.cat column
D_tall$q.type.cat = as.factor(rep(c(1:2), each = 4, times = 32))


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
# Given that there is evidence of non-normality, but NOT heteroskedasticity,
# non-parametric boostrapping and permuation
# testing will be used to estimate confidence intervals and for hypothesis checking


##############################
#### PRELIMINARY ANALYSIS ####
##############################
# analysis to determine effect of question type and 
lme.fit.prelim = lme(measure~q.type.cat+q.type+q.type.cat:q.type, random=~1|ID, data=D_tall)
anova.lme(lme.fit.prelim)


#########################
# FOLLOW UP COMPARISONS #
#########################

D_tall$q.type==0
set.seed(2018)
b = rep(0,4000) 
for(i in 1:4000){
  y = sample(D_tall$measure.2, replace=TRUE)
  lm_1 = lmer(y[D_tall$q.type==0] ~ D_tall$q.type.cat[D_tall$q.type==0] + (1|D_tall$ID[D_tall$q.type==0]), data=D_tall) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lmer(D_tall$measure.2~factor(D_tall$condition, levels=c(p,v))+(1|ID), data=D_tall)
beta_actual = fixed.effects(lm.fit)[[2]]
c(beta_actual, sum(abs(b) > beta_actual)/4000, sum(abs(b) < beta_actual)/4000,
  sum(b > beta_actual)/4000, sum(b < beta_actual)/4000)


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
  labs(x = "Test trials")