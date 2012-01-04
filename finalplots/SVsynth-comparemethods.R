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
BSMCresultsfile <- paste(resultsfolder, "BSMC/BSMC-T1000-N200000-h0.100(", sep = "")
nrunsbsmc <- 5
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
TRUEPARAMETERS <- trueparameters
NBPARAM <- dim(thetahistory)[2]
TRUESTATES <- truestates
OBSERVATIONS <- observations

# Compare the methods: possible for times 250, 500, 750 and 1000
TIME <- 500
## posterior distribution
smc2theta <- loadSMCThetaParticles(basename = SMC2Indresultsfile, nruns = nrunssmc2, time = TIME)
smc2theta$Method <- "SMC2"
bsmctheta <- loadSMCThetaParticles(basename = BSMCresultsfile, nruns = nrunsbsmc, time = TIME)
bsmctheta$Method <- "L&W"
bsmctheta$w <- 1 / 200000

PMCMCIndresultsfile <- paste(resultsfolder, "adPMCMC/adPMCMC-T", TIME, "-Iter100000-Nx500(", sep = "")
pmcmcthetas <- loadMCMCresult(basename= PMCMCIndresultsfile, nruns=1, 
                              time = 500, burnin=10000)
pmcmcthetas$Method <- "PMCMC"
thetas <- rbind(smc2theta, bsmctheta, pmcmcthetas)
thetas$Method <- factor(thetas$Method, levels=c("L&W", "SMC2", "PMCMC"))

indexparam <- 1
g <- ggplot(thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Method, nrow=1)
g <- g + xlab(explabels[indexparam])
pdf(file = paste("MethodComparTheta", indexparam, "T", TIME, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()

indexparam <- 2
g <- ggplot(thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Method, nrow=1)
g <- g + xlab(explabels[indexparam])
pdf(file = paste("MethodComparTheta", indexparam, "T", TIME, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()

indexparam <- 3
g <- ggplot(thetas, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Method, nrow=1)
g <- g + xlab(explabels[indexparam])
pdf(file = paste("MethodComparTheta", indexparam, "T", TIME, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()

indexparam <- 4
sub <- subset(thetas, Theta4 < 4)
for (m in levels(sub$Method)){
  for (t in levels(sub$Time)){
    for (r in levels(sub$Run)){
      S <- sum(sub[sub$Time == t & sub$Run == r & sub$Method == m,"w"])
      sub[sub$Time == t & sub$Run == r & sub$Method == m,"w"] <- 
        sub[sub$Time == t & sub$Run == r & sub$Method == m,"w"] / S 
    }
  }
}
g <- ggplot(sub, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Method, nrow=1)
g <- g + xlab(explabels[indexparam])
pdf(file = paste("MethodComparTheta", indexparam, "T", TIME, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()

sub <- subset(thetas, Theta5 < 0.05)
for (m in levels(sub$Method)){
  for (t in levels(sub$Time)){
    for (r in levels(sub$Run)){
      S <- sum(sub[sub$Time == t & sub$Run == r & sub$Method == m,"w"])
      sub[sub$Time == t & sub$Run == r & sub$Method == m,"w"] <- 
        sub[sub$Time == t & sub$Run == r & sub$Method == m,"w"] / S 
    }
  }
}
indexparam <- 5
g <- ggplot(sub, aes_string(x = paste("Theta", indexparam, sep = ""), weight = "w"))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Method, nrow=1)
g <- g + xlab(explabels[indexparam])
pdf(file = paste("MethodComparTheta", indexparam, "T", TIME, ".pdf", sep = ""), width = 20, height = 5)
print(g)
dev.off()
