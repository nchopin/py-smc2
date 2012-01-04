rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfile <- "/home/pierre/Dropbox/py-smc2/finalplots/LocalLevel-results/SMC2-T5000-independent-dynamicNx250-Ntheta1000(0).RData"
load(file = resultsfile)

Ntheta <- dim(thetahistory)[3]
ESSdataframe <- as.data.frame(cbind(1:length(ESS), ESS))
g <- ggplot(data = ESSdataframe, aes(x = V1, y= ESS))
g <- g + geom_line() + xlab("iterations") + ylab("ESS") + ylim(0, Ntheta)
pdf(file = "locallevel-ESS.pdf", width = 10, height = 5)
print(g)
dev.off()


g <- qplot(x = 1:T, y = cumsum(computingtimes), geom = "line",
           ylab = "computing time (square root scale)", xlab = "iteration")
g <- g + scale_y_sqrt()
pdf(file = "locallevel-computingtime.pdf", width = 10, height = 5)
print(g)
dev.off()







