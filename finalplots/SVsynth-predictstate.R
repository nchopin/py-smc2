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

## predicted state 1
smc2state1 <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "predictedstate1", 
                                  somename = c("mean", "quantile5", "quantile95"))
bsmcstate1 <- loadSomeQuantity(basename = BSMCresultsfile, nruns = nrunsbsmc,
                                  somequantity = "predictedstate1",
                                  somename = c("mean", "quantile5", "quantile95"))
smc2state1$Method <- "SMC2"
bsmcstate1$Method <- "BSMC"
state1 <- rbind(smc2state1, bsmcstate1)

start <- 10
g <- ggplot(subset(state1, Time >= start), 
            aes(x = Time, y = quantile95, colour = Method, linetype = Run))
g <- g + geom_line(alpha = 1)
g <- g + geom_line(aes(x = Time, y = TRUESTATES[start:T,1]), colour = "darkgreen", alpha = 0.2,
                   linetype = 1)
g <- g + ylab("95% quantile")
pdf(file = "predictstate.pdf", height = 5, width = 10)
print(g)
dev.off()
pdf(file = "predictstateAroundCrash.pdf", height = 5, width = 10)
print(g + xlim(400, 600))
dev.off()

sum(rep(TRUESTATES[,1], nrunsbsmc) > subset(state1, Method == "BSMC")[,"quantile95"]) / nrunsbsmc
sum(rep(TRUESTATES[,1], nrunssmc2) > subset(state1, Method == "SMC2")[,"quantile95"]) / nrunssmc2
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 1))$quantile95)
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 2))$quantile95)
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 3))$quantile95)
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 4))$quantile95)
