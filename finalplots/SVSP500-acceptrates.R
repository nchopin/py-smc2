rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

explabels <- c(expression(mu), expression(beta), expression(xi), expression(omega^2), expression(lambda))

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/SVSP500-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SVmultifactor/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
nrunssmc2 <- 3
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
NBPARAM <- dim(thetahistory)[2]
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
g <- g + xlab("Time") + xlim(0, 753) + ylab("Acceptance rates") + ylim(0, 1)
pdf(file = "SP500acceptrates.pdf")
print(g)
dev.off()



