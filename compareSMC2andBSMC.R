rm(list = ls())
gc()
graphics.off()
library(ggplot2)
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/SVonefactor/synthetic/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2-T100-independent-dynamicNx100-Ntheta1000(", sep = "")
nrunssmc2 <- 3
BSMCresultsfile <- paste(resultsfolder, "BSMC-T100-N20000-h0.100(", sep = "")
nrunsbsmc <- 3
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
TRUEPARAMETERS <- trueparameters
NBPARAM <- dim(thetahistory)[2]
TRUESTATES <- truestates
OBSERVATIONS <- observations

## predicted squared observations
smc2predsqobs <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "predictedsquaredobs",
                                  somename = c("mean", "quantile5", "quantile95"))
bsmcpredsqobs <- loadSomeQuantity(basename = BSMCresultsfile, nruns = nrunsbsmc,
                                  somequantity = "predictedsquaredobs",
                                  somename = c("mean", "quantile5", "quantile95"))

smc2predsqobs$Method <- "SMC2"
bsmcpredsqobs$Method <- "BSMC"
predsqobs <- rbind(smc2predsqobs, bsmcpredsqobs)

g <- ggplot(subset(predsqobs, Time > 10), aes(x = Time, y = mean, colour = Run, linetype = "mean"))
g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + geom_line(aes(y = quantile5, linetype = "quantile"))
g <- g + geom_line(aes(y = quantile95, linetype = "quantile"))
g <- g + ylab("squared observations") + xlab("time") + scale_linetype_discrete(name = "Quantity")
#print(g)
## log evidences
smc2evid <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "evidences", somename = "LogEvid",
                                  functiontoapply = cumsum)
bsmcevid <- loadSomeQuantity(basename = BSMCresultsfile, nruns = nrunsbsmc,
                                  somequantity = "evidences", somename = "LogEvid", 
                                  functiontoapply = cumsum)
smc2evid$Method <- "SMC2"
bsmcevid$Method <- "BSMC"
evid <- rbind(smc2evid, bsmcevid)
g <- ggplot(evid, aes(x = Time, y = LogEvid, colour = Run))
g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + ylab("evid") + xlab("time")
#print(g)
## computing time
smc2CT <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "computingtimes", somename = "CT",
                                  functiontoapply = cumsum)
bsmcCT <- loadSomeQuantity(basename = BSMCresultsfile, nruns = nrunsbsmc,
                                  somequantity = "computingtimes", somename = "CT", 
                                  functiontoapply = cumsum)
smc2CT$Method <- "SMC2"
bsmcCT$Method <- "BSMC"
CT <- rbind(smc2CT, bsmcCT)
g <- ggplot(CT, aes(x = Time, y = CT, colour = Run))
g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + ylab("computing time") + xlab("time")
#print(g)

TIME <- c(100)
## posterior distribution
smc2theta <- loadSMCThetaParticles(basename = SMC2Indresultsfile, nruns = nrunssmc2, time = TIME)
smc2theta$Method <- "SMC2"
bsmctheta <- loadSMCThetaParticles(basename = BSMCresultsfile, nruns = nrunsbsmc, time = TIME)
bsmctheta$Method <- "BSMC"
thetas <- rbind(smc2theta, bsmctheta)

indexparam <- 3
g <- ggplot(thetas, aes_string(x = paste("Theta", indexparam, sep = "")))
g <- g + geom_density(aes(colour = Run), alpha = 0.2)
g <- g + facet_wrap(~ Method)
g <- g + geom_vline(xintercept = TRUEPARAMETERS[indexparam], linetype = 3)
print(g)


