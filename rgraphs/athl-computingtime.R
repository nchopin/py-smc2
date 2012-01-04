rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/Athletics-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2-T35-independent-dynamicNx250-Ntheta1000(", sep = "")
nrunssmc2 <- 10

## computing time
smc2CT <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "computingtimes", somename = "CT",
                                  functiontoapply = cumsum)
g <- ggplot(smc2CT, aes(x = Time, y = CT, colour = Run))
g <- g + geom_line()
g <- g + ylab("computing time") + xlab("time")
pdf(file = "Athl-ComputingTime.pdf")
print(g)
dev.off()
