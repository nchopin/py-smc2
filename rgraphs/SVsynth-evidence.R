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
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
TRUEPARAMETERS <- trueparameters
NBPARAM <- dim(thetahistory)[2]
TRUESTATES <- truestates
OBSERVATIONS <- observations


## log evidences
smc2evid <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "evidences", somename = "LogEvid",
                                  functiontoapply = function(x) cumsum(log(x)))
bsmcevid <- loadSomeQuantity(basename = BSMCresultsfile, nruns = nrunsbsmc,
                                  somequantity = "evidences", somename = "LogEvid", 
                                  functiontoapply = function(x) cumsum(log(x)))
smc2evid$Method <- "SMC2"
bsmcevid$Method <- "L&W"
evid <- rbind(smc2evid, bsmcevid)


#head(evid)
meanevid <- (aggregate(x = evid$LogEvid, by = list(evid$Time, evid$Method), FUN = mean))
names(meanevid) <- c("Time", "Method", "mean")
#g <- ggplot(subset(meanevid, Time < 410 & Time > 400), 
#            aes(x = Time, y = mean, colour = Method))
#g <- g + geom_line() + ylab("log evidence (mean across runs)")
#pdf(file = "LogEvidenceAroundCrash.pdf")
#print(g)
#dev.off()
#
diffevid <- aggregate(x = meanevid$mean, by = list(meanevid$Time), FUN = function(x) x[1] - x[2])
#head(diffevid)
names(diffevid) <- c("Time", "diff")
#g <- ggplot(subset(diffevid, Time < T & Time > 0), 
#            aes(x = Time, y = diff))
#g <- g + geom_line() + ylab("log evidence: difference between methods")
#pdf(file = "LogEvidenceDiffBetweenMethods.pdf")
#print(g)
#dev.off()

base <- subset(meanevid, Method == "SMC2")$mean
head(evid)
evid$LogEvidMinusBase <- evid$LogEvid - rep(base, nrunsbsmc + nrunssmc2)
g <- ggplot(subset(evid, Time < T & Time > 10), 
            aes(x = Time, y = LogEvidMinusBase, colour = Run))
g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + ylab("log evidence minus reference") + xlab("Time")
pdf(file = "LogEvidenceVariationsAcrossRuns.pdf", width = 15, height = 5)
print(g)
dev.off()



