rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

explabels <- c(expression(mu), expression(beta), expression(xi), expression(omega^2), expression(lambda))

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/SVsynthetic-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2/SMC2-T1000-independent-dynamicNx100-Ntheta1000(", sep = "")
nrunssmc2 <- 5
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
TRUEPARAMETERS <- trueparameters
NBPARAM <- dim(thetahistory)[2]
TRUESTATES <- truestates
OBSERVATIONS <- observations


## posterior distribution
smc2theta250 <- loadSMCThetaParticles(basename = SMC2Indresultsfile, 
                                   nruns = nrunssmc2, time = 250)
smc2theta500 <- loadSMCThetaParticles(basename = SMC2Indresultsfile, 
                                   nruns = nrunssmc2, time = 500)
smc2theta1000 <- loadSMCThetaParticles(basename = SMC2Indresultsfile, 
                                   nruns = nrunssmc2, time = 1000)
smc2thetas <- rbind(smc2theta250, smc2theta500, smc2theta1000)
head(smc2thetas)

indexparam <- 1
g <- ggplot(smc2thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Time, nrow=1)
g <- g + geom_vline(xintercept = TRUEPARAMETERS[indexparam], linetype = 3)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dnorm(x, sd = 1.41421)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
pdf(file = paste("SMC2theta", indexparam, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()
indexparam <- 2
g <- ggplot(smc2thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Time, nrow=1)
g <- g + geom_vline(xintercept = TRUEPARAMETERS[indexparam], linetype = 3)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dnorm(x, sd = 1.41421)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
pdf(file = paste("SMC2theta", indexparam, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()

indexparam <- 3
g <- ggplot(smc2thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Time, nrow=1)
g <- g + geom_vline(xintercept = TRUEPARAMETERS[indexparam], linetype = 3)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
pdf(file = paste("SMC2theta", indexparam, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()

indexparam <- 4
sub <- subset(smc2thetas, Theta4 < 4)
for (t in levels(sub$Time)){
  for (r in levels(sub$Run)){
    S <- sum(sub[sub$Time == t & sub$Run == r,"w"])
    sub[sub$Time == t & sub$Run == r,"w"] <- sub[sub$Time == t & sub$Run == r,"w"] / S 
  }
}
g <- ggplot(sub, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Time, nrow=1)
g <- g + geom_vline(xintercept = TRUEPARAMETERS[indexparam], linetype = 3)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
pdf(file = paste("SMC2theta", indexparam, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()
####### need to normalize the weights so that the density sum to 1

sub <- subset(smc2thetas, Theta5 < 0.05)
for (t in levels(sub$Time)){
  for (r in levels(sub$Run)){
    S <- sum(sub[sub$Time == t & sub$Run == r,"w"])
    sub[sub$Time == t & sub$Run == r,"w"] <- sub[sub$Time == t & sub$Run == r,"w"] / S 
  }
}
indexparam <- 5
g <- ggplot(sub, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Time, nrow=1)
g <- g + geom_vline(xintercept = TRUEPARAMETERS[indexparam], linetype = 3)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 1)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
pdf(file = paste("SMC2theta", indexparam, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()


