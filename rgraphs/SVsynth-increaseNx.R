rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/SVsynthetic-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2/SMC2-T1000-independent-dynamicNx100-Ntheta1000(", sep = "")
nrunssmc2 <- 5
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
TRUEPARAMETERS <- trueparameters
NBPARAM <- dim(thetahistory)[2]
TRUESTATES <- truestates
OBSERVATIONS <- observations

library(foreach)
dfNx <- foreach(run = 1:nrunssmc2, .combine = rbind) %do% {
    load(file = paste(SMC2Indresultsfile, run - 1, ").RData", sep = ""))
    Nxlist <- c(Nxlist, Nxlist[length(Nxlist)])
    increaseindices <- c(increaseindices, T)
    quantityDF <- as.data.frame(cbind(Nxlist, increaseindices, run))
#     names(quantityDF) <- c("acceptratios", "indices", "Run")
    quantityDF
}

dfNx$run <- as.factor(dfNx$run)

g <- ggplot(dfNx, aes(x = increaseindices, y = Nxlist, colour = run))
g <- g + geom_step() + geom_point()
g <- g + xlab("Time") + xlim(0, 1000) + ylab(expression(N[x]))
pdf(file = "SMC2increaseNx.pdf")
print(g)
dev.off()


