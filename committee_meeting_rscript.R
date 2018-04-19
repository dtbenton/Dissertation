# Create a function to plot significance stars and lines
sig_func = function(y,offset,a,b){
  y <- y
  # set an offset for tick lengths
  offset <- offset
  # draw first horizontal line
  lines(x[c(a,b)],c(y, y))
  # draw ticks
  lines(x[c(a,a)],c(y, y-offset))
  lines(x[c(b,b)],c(y, y-offset))
  # draw asterics
  text(x[a]+((x[b]-x[a])/b),y+offset,"*")
}

### EXP 1 - ADULTS (n = 60) ###
options(scipen=999)
D = read.csv(file.choose(), header=TRUE)

# install relevant packages
install.packages('lazyeval')
install.packages('ggplot2', dep=TRUE)


# reshape and reorder data
D_tall = reshape(D, varying = 2:9, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:512, direction = "long")
D_tall = D_tall[order(D_tall$ID),]
D_tall$condition = as.factor(D_tall$condition)
head(D_tall)


## barplot for perceptual ratings ##
perceptual.bars = c((100-mean(D_tall$measure[D_tall$condition==1])), 
                    (100-mean(D_tall$measure[D_tall$condition==7])), 
                    (100-mean(D_tall$measure[D_tall$condition==3])), 
                    (100-mean(D_tall$measure[D_tall$condition==5]))) 
x = barplot(perceptual.bars, names.arg = c("GBGR", "GRGB", "GBRG", "GRBG"), density = rep(10,4), ylim = c(0,50), main = "Exp 1 - Perceptual Ratings", 
            ylab = "Perceptual Ratings")


# PLOT SIGNIFICANCE STARS AND LINES
sig_func(42, 0.9, 1,3)
sig_func(44, 0.9, 1,4)
sig_func(38, 0.9, 2,3)
sig_func(40, 0.9, 2,4)




## barplot for causal ratings ##
causal.bars = c((100-mean(D_tall$measure[D_tall$condition==2])), 
                    (100-mean(D_tall$measure[D_tall$condition==8])), 
                    (100-mean(D_tall$measure[D_tall$condition==4])), 
                    (100-mean(D_tall$measure[D_tall$condition==6]))) 
x = barplot(causal.bars, names.arg = c("GBGR", "GRGB", "GBRG", "GRBG"), density = rep(10,4), ylim = c(0,30), main = "Exp 1 - Causal Ratings",
            ylab = "Causal Ratings")


# 2 v 8 - sig
wilcox.test(D_tall$measure[D_tall$condition==2],D_tall$measure[D_tall$condition==8], paired=TRUE)
# 2 v 4
wilcox.test(D_tall$measure[D_tall$condition==2],D_tall$measure[D_tall$condition==4], paired=TRUE)
# 2 v 6
wilcox.test(D_tall$measure[D_tall$condition==2],D_tall$measure[D_tall$condition==6], paired=TRUE)
# 8 v 4 - sig
wilcox.test(D_tall$measure[D_tall$condition==4],D_tall$measure[D_tall$condition==8], paired=TRUE)
# 8 v 6 - sig
wilcox.test(D_tall$measure[D_tall$condition==6],D_tall$measure[D_tall$condition==8], paired=TRUE)


# PLOT SIGNIFICANCE STARS AND LINES
sig_func(22, 0.9, 1,2)
sig_func(26, 0.9, 2,3)
sig_func(28, 0.9, 2,4)





###############################
### EXP 2 - ADULTS (n = 35) ###
###############################

options(scipen=999)
D = read.csv(file.choose(), header=TRUE)

# reshape and reorder data
D_tall = reshape(D, varying = 2:9, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:264, direction = "long")
D_tall = D_tall[order(D_tall$ID),]
D_tall$condition = as.factor(D_tall$condition)
head(D_tall)



## barplot for perceptual ratings ##
perceptual.bars = c((100-mean(D_tall$measure[D_tall$condition==1])), 
                    (100-mean(D_tall$measure[D_tall$condition==5])), 
                    (100-mean(D_tall$measure[D_tall$condition==3])), 
                    (100-mean(D_tall$measure[D_tall$condition==7]))) 
x = barplot(perceptual.bars, names.arg = c("GBGR", "GBgapGR", "GBRG", "GBgapRG"), density = rep(10,4), ylim = c(0,60), main = "Exp 2 - Perceptual Ratings", 
            ylab = "Perceptual Ratings")


# paired t-tests
# 1 v 2 - sig
wilcox.test(D_tall$measure[D_tall$condition==1],D_tall$measure[D_tall$condition==5], paired=TRUE)

# 1 v 3
wilcox.test(D_tall$measure[D_tall$condition==1],D_tall$measure[D_tall$condition==3], paired=TRUE)

# 1 v 4 - sig
wilcox.test(D_tall$measure[D_tall$condition==1],D_tall$measure[D_tall$condition==7], paired=TRUE)

# 2 v 4
wilcox.test(D_tall$measure[D_tall$condition==5],D_tall$measure[D_tall$condition==7], paired=TRUE)

# 2 v 3 - sig
wilcox.test(D_tall$measure[D_tall$condition==5],D_tall$measure[D_tall$condition==3], paired=TRUE)

# 3 v 4
wilcox.test(D_tall$measure[D_tall$condition==3],D_tall$measure[D_tall$condition==7], paired=TRUE)


# PLOT SIGNIFICANCE STARS AND LINES
sig_func(51, 0.9, 1,2)
sig_func(54, 0.9, 1,4)
sig_func(51, 0.9, 3,4)



## barplot for causal ratings ##
causal.bars = c((100-mean(D_tall$measure[D_tall$condition==2])), 
                    (100-mean(D_tall$measure[D_tall$condition==6])), 
                    (100-mean(D_tall$measure[D_tall$condition==4])), 
                    (100-mean(D_tall$measure[D_tall$condition==8]))) 
x = barplot(causal.bars, names.arg = c("GBGR", "GBgapGR", "GBRG", "GBgapRG"), density = rep(10,4), ylim = c(0,60), main = "Exp 2 - Causal Ratings", 
            ylab = "Causal Ratings")


# paired t-tests
# 1 v 2 - sig
wilcox.test(D_tall$measure[D_tall$condition==2],D_tall$measure[D_tall$condition==6], paired=TRUE)

# 1 v 3
wilcox.test(D_tall$measure[D_tall$condition==2],D_tall$measure[D_tall$condition==4], paired=TRUE)

# 1 v 4 - sig
wilcox.test(D_tall$measure[D_tall$condition==2],D_tall$measure[D_tall$condition==8], paired=TRUE)

# 2 v 4
wilcox.test(D_tall$measure[D_tall$condition==6],D_tall$measure[D_tall$condition==8], paired=TRUE)

# 2 v 3 - sig
wilcox.test(D_tall$measure[D_tall$condition==6],D_tall$measure[D_tall$condition==4], paired=TRUE)

# 3 v 4
wilcox.test(D_tall$measure[D_tall$condition==3],D_tall$measure[D_tall$condition==7], paired=TRUE)


# PLOT SIGNIFICANCE STARS AND LINES
sig_func(51, 0.9, 1,2)
sig_func(54, 0.9, 1,4)
sig_func(51, 0.9, 3,4)






###############################
### EXP 1 - Infants (n = 35) ###
###############################
D = read.csv(file.choose(), header=TRUE)


# reshape and reorder data
D_tall = reshape(D, varying = 2:5, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:65, direction = "long")
D_tall = D_tall[order(D_tall$ID),]
D_tall$condition = as.factor(D_tall$condition)


causal.bars = c(mean(D_tall$measure[D_tall$condition==1]), 
                mean(D_tall$measure[D_tall$condition==2]), 
                mean(D_tall$measure[D_tall$condition==3]), 
                mean(D_tall$measure[D_tall$condition==4])) 
x = barplot(causal.bars, names.arg = c("GBGR", "GRGB", "GBRG", "GRBG"), density = rep(10,4), ylim = c(0,15), main = "Exp 1 - Infants",
            ylab = "Looking Time (s)")


# paired t-tests
wilcox_func = function(a,b){
  calc = wilcox.test(D_tall$measure[D_tall$condition==a],D_tall$measure[D_tall$condition==b], paired=TRUE)
  calc$p.value
}


wilcox_func(1,2)
wilcox_func(1,3)
wilcox_func(1,4)
wilcox_func(2,3)
wilcox_func(2,4)
wilcox_func(3,4)