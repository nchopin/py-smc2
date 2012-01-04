rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/SVSP500-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SVmultifactor/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
nrunssmc2 <- 3
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
g <- g + xlab("Time") + xlim(0, 753) + ylab(expression(N[x]))
pdf(file = "SP500increaseNx.pdf")
print(g)
dev.off()



