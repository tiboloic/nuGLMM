# LT 30/10

# a temporary list of tests of nuglmm
# FIXME: implement proper tests to be performed after building


# simple intercept and slope with random effect
x = rnorm(100)
grp = gl(20,5)
X = cbind(1, x)
beta = c(1,1)
b = rnorm(20)
Z = apply(diag(20), 1, function(l) t(replicate(5,l)))
y = rpois(100, exp(X %*% beta + Z %*% b))
Y = replicate(5, rpois(100, exp(X %*% beta + Z %*% b)))

require(glmmTMB)
#m = glmmTMB(y~x + (1|grp), family=poisson())
#m2 = nuglmm(Y~x + (1|grp), family=poisson())

# fist inference test: a random effect 
m1 = nuglmm(Y~x, family=poisson())
m2 = nuglmm(Y~x + (1|grp), family="poisson")

# test offset
m1 = nuglmm(Y~x, family=poisson(log), offset=1)

# test matrix offset
m1 = nuglmm(Y~x, family=poisson(), offset=matrix(1, 100, 5))

# test dummy familly
m1 = nuglmm(Y~x, family=tiboloid)

# variable lengths in fomula
m1 = nuglmm(Y~rep(x,2), family=poisson())
            
# try presence absence data
sum(Y==0)
Y2 = apply(Y, c(1,2), function(y) ifelse(y>0, 1, 0))
m1 = nuglmm(Y2~x, family=binomial(cloglog))
m2 = nuglmm(Y2~x + (1|grp), family=binomial(cloglog))
anova(m1, m2)
system.time(print(nuglmm.anova(m1,m2, nboot=10)))
#m = glmmTMB(y~x + (1|grp), family=poisson())
#m2 = nuglmm(Y~x + (1|grp), family=poisson())

# fist inference test: a random effect 
m1 = nuglmm(Y~x, family=poisson())
m2 = nuglmm(Y~x + (1|grp), family="poisson")
anova(m1,m2)
system.time(print(nuglmm.anova(m1,m2, nboot=500)))




###############################################################################################################
# Tasmania

library(mvabund)
data(Tasmania)

# first test distribution
m1 = nuglmm(abund ~ treatment * block, data=Tasmania, family="poisson")
m2 = nuglmm(abund ~ treatment * block, data=Tasmania, family="nbinom2")
anova(m1, m2)
# FIXME BUG
#system.time(print(nuglmm.anova(m1,m2, nboot=50)))

# test for block
m1 = nuglmm(abund ~ treatment, data=Tasmania, family="nbinom2")
m2 = nuglmm(abund ~ treatment + block, data=Tasmania, family="nbinom2")
anova(m1, m2)
# p = 0.01 , 20s 
system.time(print(nuglmm.anova(m1,m2, nboot=50)))

# test for treat
m1 = nuglmm(abund ~ block, data=Tasmania, family="nbinom2")
m2 = nuglmm(abund ~ treatment + block, data=Tasmania, family="nbinom2")
anova(m1, m2)
# p = 0.01 , 25s 
nuglmm.anova(m1,m2)

# test for fixed effect interaction between treatment and block
m1 = nuglmm(abund ~ treatment + block, data=Tasmania, family="nbinom2")
m2 = nuglmm(abund ~ treatment * block, data=Tasmania, family="nbinom2")

anova(m1,m2)
nuglmm.anova(m1,m2)
# p = 0.08, 25s
system.time(print(nuglmm.anova(m1,m2, nboot=500)))
# p = 0.10 with 500 simul

# compare with mvabund
tas = as.mvabund(Tasmania$abund)
m1. = manyglm(tas ~ treatment * block, data=Tasmania, family="negative.binomial")
summary(m1., nBoot=100)
anova(m1., nBoot=100)

# we want block to be a random effect, test for treatment:block two formulations
# formulation 1

m1 = nuglmm(abund ~ treatment + (1|block), data=Tasmania, family="nbinom2")
m2 = nuglmm(abund ~ treatment + (1|block) + (1|treatment:block), data=Tasmania, family="nbinom2")
anova(m1,m2)
# 10 minutes for 100 p = 0.11
system.time(print(nuglmm.anova(m1,m2, nboot=100)))


# formulation 2

m1 = nuglmm(abund ~ treatment + (1|block), data=Tasmania, family="nbinom2")
m2 = nuglmm(abund ~ treatment + (treatment|block), data=Tasmania, family="nbinom2")
anova(m1,m2)
# 10 minutes for 100 p = 0.04 !
system.time(print(nuglmm.anova(m1,m2, nboot=50)))

############################################################################################
# spider
############################################################################################
data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

#To fit a log-linear model assuming counts are poisson:
glm.spid <- manyglm(spiddat~X, family="poisson")
glm.spid 

summary(glm.spid, resamp="residual")

spiddat = as.data.frame(spider$x)
spiddat$abund = as.matrix(spider$abund)
m1 = nuglmm(abund~soil.dry + fallen.leaves + moss + herb.layer + reflection, data = spiddat, family="poisson")
m2 = nuglmm(abund~soil.dry + bare.sand + fallen.leaves + moss + herb.layer + reflection, data = spiddat, family="poisson")
# FIXME
anova(m1, m2)


