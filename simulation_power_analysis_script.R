knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE)

# define the parameters
mu = c(80, 76, 77, 83) # true effects (in this case, a double dissociation)
sigma = 20  # population standard deviation
rho = 0.75 # correlation between repeated measures
nsubs = 16 # how many subjects?
nsims = 5000 # how many simulation replicates?

# create 2 factors representing the 2 independent variables
cond = data.frame(
  X1 = rep(factor(letters[1:4]), nsubs * 2))

# create a subjects factor
subject = factor(sort(rep(1:nsubs, 4)))

# combine above into the design matrix
dm = data.frame(subject, cond)


# create k x k matrix populated with sigma
sigma.mat <- rep(sigma, 4)
S <- matrix(sigma.mat, ncol=length(sigma.mat), nrow=length(sigma.mat))

# compute covariance between measures
Sigma <- t(S) * S * rho  

# put the variances on the diagonal 
diag(Sigma) <- sigma^2  






# stack 'nsims' individual data frames into one large data frame
df = dm[rep(seq_len(nrow(dm)), nsims), ]

# add an index column to track the simulation run
df$simID = sort(rep(seq_len(nsims), nrow(dm)))

# sample the observed data from a multivariate normal distribution
# using MASS::mvrnorm with the parameters mu and Sigma created earlier
# and bind to the existing df

require(MASS)
make.y = expression(as.vector(t(mvrnorm(nsubs, mu, Sigma))))
df$y = as.vector(replicate(nsims, eval(make.y)))             

# use do(), the general purpose complement to the specialized data 
# manipulation functions available in dplyr, to run the ANOVA on
# each section of the grouped data frame created by group_by

require(dplyr)
require(car)
require(broom)

mods <- df %>% 
  group_by(simID) %>% 
  do(model = aov(y ~ X1 + Error(subject / (X1)), qr=FALSE, data = .)) 


# extract p-values for each effect and store in a data frame
p = data.frame(
  mods %>% do(as.data.frame(tidy(.$model[[3]])$p.value[1])))
colnames(p) = c('X1')


power = apply(as.matrix(p), 2, 
              function(x) round(mean(ifelse(x < .05, 1, 0) * 100),0))