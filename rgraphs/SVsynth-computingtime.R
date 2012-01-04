rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/SVsynthetic-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2/SMC2-T1000-independent-dynamicNx100-Ntheta1000(", sep = "")
nrunssmc2 <- 5
BSMCresultsfile <- paste(resultsfolder, "BSMC/BSMC-T1000-N200000-h0.100(", sep = "")
nrunsbsmc <- 5

## computing time
smc2CT <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "computingtimes", somename = "CT",
                                  functiontoapply = cumsum)
bsmcCT <- loadSomeQuantity(basename = BSMCresultsfile, nruns = nrunsbsmc,
                                  somequantity = "computingtimes", somename = "CT", 
                                  functiontoapply = cumsum)
smc2CT$Method <- "SMC2"
bsmcCT$Method <- "L&W"
CT <- rbind(smc2CT, bsmcCT)
g <- ggplot(CT, aes(x = Time, y = CT, colour = Run))
g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + ylab("computing time") + xlab("time")
pdf(file = "ComputingTime.pdf", width = 10, height = 5)
print(g)
dev.off()
