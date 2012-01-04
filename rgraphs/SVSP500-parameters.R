rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/SVSP500-results/"
singlefactorfile <- paste(resultsfolder, "SVonefactor/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
twofactorfile <- paste(resultsfolder, "SVfixedrho/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
fullfile <- paste(resultsfolder, "SVmultifactor/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
nrunssmc2 <- 3

onefactortheta <- loadSMCThetaParticles(basename = singlefactorfile, 
                                   nruns = nrunssmc2, time = 753)
onefactortheta$Model <- "SV one factor"
twofactortheta <- loadSMCThetaParticles(basename = twofactorfile, 
                                   nruns = nrunssmc2, time = 753)
twofactortheta$Model <- "SV multifactor"
fulltheta <- loadSMCThetaParticles(basename = fullfile, 
                                   nruns = nrunssmc2, time = 753)
fulltheta$Model <- "SV full"

explabels <- c(expression(mu), expression(beta), expression(xi), expression(omega^2), expression(lambda[1]))
columnstokeep <- c(paste("Theta", 1:5, sep = ""), "Run", "Time", "w", "Model")
theta <- onefactortheta[, names(onefactortheta) %in% columnstokeep]
theta <- rbind(theta, twofactortheta[, names(twofactortheta) %in% columnstokeep])
theta <- rbind(theta, fulltheta[, names(fulltheta) %in% columnstokeep])

indexparam <- 1
g <- ggplot(theta, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(fill = Model), alpha = 0.8)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dnorm(x, sd = 1.41421)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + scale_fill_discrete(name = "Model:", h = c(0, 0, 0), c = c(0, 0, 0), l = c(0, 50, 100))
g <- g + opts(legend.position = "none")
pdf(file = paste("SP500Theta", indexparam, ".pdf", sep = ""))
print(g)
dev.off()

indexparam <- 2
g <- ggplot(theta, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(fill = Model), alpha = 0.8)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dnorm(x, sd = 1.41421)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + scale_fill_discrete(name = "Model:", h = c(0, 0, 0), c = c(0, 0, 0), l = c(0, 50, 100))
g <- g + opts(legend.position = "none")
pdf(file = paste("SP500Theta", indexparam, ".pdf", sep = ""))
print(g)
dev.off()

indexparam <- 3
g <- ggplot(theta, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(fill = Model), alpha = 0.8)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + scale_fill_discrete(name = "Model:", h = c(0, 0, 0), c = c(0, 0, 0), l = c(0, 50, 100))
g <- g + opts(legend.position = "none")
pdf(file = paste("SP500Theta", indexparam, ".pdf", sep = ""))
print(g)
dev.off()

indexparam <- 4
g <- ggplot(theta, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(fill = Model), alpha = 0.8)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + scale_fill_discrete(name = "Model:", h = c(0, 0, 0), c = c(0, 0, 0), l = c(0, 50, 100))
g <- g + opts(legend.position = "none")
pdf(file = paste("SP500Theta", indexparam, ".pdf", sep = ""))
print(g)
dev.off()

indexparam <- 5
g <- ggplot(theta, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(fill = Model), alpha = 0.8)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 1)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + scale_fill_discrete(name = "Model:", h = c(0, 0, 0), c = c(0, 0, 0), l = c(0, 50, 100))
g <- g + opts(legend.position = "none")
pdf(file = paste("SP500Theta", indexparam, ".pdf", sep = ""))
print(g)
dev.off()

columnstokeep <- c(paste("Theta", 6:7, sep = ""), "Run", "Time", "w", "Model")
explabels <- c(expression(lambda[2] - lambda[1]), expression(w[1]))
theta <- twofactortheta[, names(twofactortheta) %in% columnstokeep]
theta <- rbind(theta, fulltheta[, names(fulltheta) %in% columnstokeep])

indexparam <- 1
g <- ggplot(theta, aes_string(x = paste("Theta", 5 + indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(fill = Model), alpha = 0.8)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) dexp(x, rate = 0.5)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + scale_fill_discrete(name = "Model:", h = c(0, 0, 0), c = c(0, 0, 0), l = c(0, 50, 100))
g <- g + opts(legend.position = "none")
pdf(file = paste("SP500Theta", 5 + indexparam, ".pdf", sep = ""))
print(g)
dev.off()

indexparam <- 2
g <- ggplot(theta, aes_string(x = paste("Theta", 5 + indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(fill = Model), alpha = 0.8)
g <- g + xlab(explabels[indexparam])
priorfunction <- function(x) 1
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + scale_fill_discrete(name = "Model:", h = c(0, 0, 0), c = c(0, 0, 0), l = c(0, 50, 100))
g <- g + opts(legend.position = "none")
pdf(file = paste("SP500Theta", 5 + indexparam, ".pdf", sep = ""))
print(g)
dev.off()


