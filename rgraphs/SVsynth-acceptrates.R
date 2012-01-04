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

library(foreach)
dfacc <- foreach(run = 1:nrunssmc2, .combine = rbind) %do% {
    load(file = paste(SMC2Indresultsfile, run - 1, ").RData", sep = ""))
    quantityDF <- as.data.frame(cbind(acceptratios, resamplingindices, run))
#     names(quantityDF) <- c("acceptratios", "indices", "Run")
    quantityDF
}

dfacc$run <- as.factor(dfacc$run)
dfacc
g <- ggplot(dfacc, aes(x = resamplingindices, y = acceptratios, colour = run))
g <- g + geom_line() + geom_point()
g <- g + xlab("Time") + xlim(0, 1000) + ylab("Acceptance rates") + ylim(0, 1)
pdf(file = "SMC2acceptrates.pdf")
print(g)
dev.off()



