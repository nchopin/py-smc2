rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

explabels <- c(expression(nu), expression(xi), expression(sigma))

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/Athletics-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2-T35-independent-dynamicNx250-Ntheta1000(", sep = "")
nrunssmc2 <- 10
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
NBPARAM <- dim(thetahistory)[2]
OBSERVATIONS <- observations


## posterior distribution
smc2thetas <- loadSMCThetaParticles(basename = SMC2Indresultsfile, 
                                   nruns = nrunssmc2, time = 35)
head(smc2thetas)

indexparam <- 1
g <- ggplot(smc2thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run, fill = Run), alpha = 0.2)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
pdf(file = "Athl-Theta1.pdf")
print(g)
dev.off()

indexparam <- 2
g <- ggplot(smc2thetas, aes_string(x = paste("-Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run, fill = Run), alpha = 0.2)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(-x, rate = 0.50000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
pdf(file = "Athl-Theta2.pdf")
print(g)
dev.off()

indexparam <- 3
g <- ggplot(smc2thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run, fill = Run), alpha = 0.2)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
pdf(file = "Athl-Theta3.pdf")
print(g)
dev.off()


