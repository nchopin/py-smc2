rm(list = ls())
gc()
graphics.off()
library(ggplot2)
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/SVonefactor/synthetic/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2-T1000-independent-dynamicNx100-Ntheta1000(", sep = "")
nrunssmc2 <- 5
BSMCresultsfile <- paste(resultsfolder, "BSMC-T1000-N150000-h0.100(", sep = "")
nrunsbsmc <- 5
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
sqOBS <- OBSERVATIONS**2
start <- 10
g <- ggplot(subset(predsqobs, Time >= start), aes(x = Time, y = mean, colour = Run, linetype = "mean"))
g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + geom_line(aes(x = Time, y = sqOBS[start:T]), colour = "red", alpha = 0.4)
g <- g + geom_line(aes(y = quantile5, linetype = "quantile"))
g <- g + geom_line(aes(y = quantile95, linetype = "quantile"))
g <- g + ylab("squared observations") + xlab("time") + scale_linetype_discrete(name = "Quantity")
# print(g)
## other representation
start <- 10
g <- ggplot(subset(predsqobs, Time >= start), 
            aes(x = Time, y = quantile95, colour = Method, linetype = Run))
g <- g + geom_line(alpha = 1)
g <- g + geom_line(aes(x = Time, y = sqOBS[start:T]), colour = "darkgreen", alpha = 0.2,
                   linetype = 1)
# pdf(file = "SquaredObsPrediction.pdf", width = 10, height = 5)
# print(g)
# dev.off()

sum(rep(sqOBS, nrunsbsmc) > subset(predsqobs, Method == "BSMC")[,"quantile95"]) / nrunsbsmc
sum(rep(sqOBS, nrunssmc2) > subset(predsqobs, Method == "SMC2")[,"quantile95"]) / nrunssmc2

## log evidences
smc2evid <- loadSomeQuantity(basename = SMC2Indresultsfile, nruns = nrunssmc2,
                                  somequantity = "evidences", somename = "LogEvid",
                                  functiontoapply = function(x) cumsum(log(x)))
bsmcevid <- loadSomeQuantity(basename = BSMCresultsfile, nruns = nrunsbsmc,
                                  somequantity = "evidences", somename = "LogEvid", 
                                  functiontoapply = function(x) cumsum(log(x)))
smc2evid$Method <- "SMC2"
bsmcevid$Method <- "BSMC"
evid <- rbind(smc2evid, bsmcevid)
g <- ggplot(subset(evid, Time < 100 & Time > 80), aes(x = Time, y = LogEvid, colour = Run, linetype = Method))
# g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + ylab("cumulated log evidence") + xlab("time")
print(g)

# pdf(file = "CumulatedLogEvidenceAroundShock.pdf", width = 10, height = 5)
# print(g)
# print(g + xlim(400, 410) + ylim(-450, -420))
head(evid)
pdf(file = "EvidenceComparisonSMC2vsBSMC.pdf")
meanevid <- (aggregate(x = evid$LogEvid, by = list(evid$Time, evid$Method), FUN = mean))
names(meanevid) <- c("Time", "Method", "mean")
g <- ggplot(subset(meanevid, Time < 410 & Time > 400), 
            aes(x = Time, y = mean, colour = Method))
g <- g + geom_line() + ylab("log evidence (mean across runs)")
print(g)

diffevid <- aggregate(x = meanevid$mean, by = list(meanevid$Time), FUN = function(x) x[1] - x[2])
head(diffevid)
names(diffevid) <- c("Time", "diff")
g <- ggplot(subset(diffevid, Time < T & Time > 0), 
            aes(x = Time, y = diff))
g <- g + geom_line() + ylab("difference between means [across runs]")
print(g)

base <- subset(meanevid, Method == "SMC2")$mean
head(evid)
evid$LogEvidMinusBase <- evid$LogEvid - rep(base, nrunsbsmc + nrunssmc2)
g <- ggplot(subset(evid, Time < T & Time > 10), 
            aes(x = Time, y = LogEvidMinusBase, colour = Method, linetype = Run))
# g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + ylab("log evidence minus base (=mean of SMC2 runs)") + xlab("time")
print(g)
dev.off()

# print(g + xlim(950, 1000) + ylim(-1000, -900))
# dev.off()
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
# pdf(file = "ComputingTime.pdf", width = 10, height = 5)
# print(g)
# dev.off()

TIME <- c(500)
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
# pdf(file = "ThirdParameterAtT500.pdf", width = 10, height = 5)
# print(g)
# dev.off()

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
start <- 1
g <- ggplot(subset(state1, Time >= start), aes(x = Time, y = quantile5, colour = Run, 
                                               linetype = "quantile"))
g <- g + facet_wrap(~ Method)
g <- g + geom_line()
g <- g + geom_line(aes(x = Time, y = TRUESTATES[start:T,1]), colour = "red", alpha = 0.4)
# g <- g + geom_line(aes(y = quantile5, linetype = "quantile"))
g <- g + geom_line(aes(y = quantile95, linetype = "quantile"))
g <- g + ylab("state 1") + xlab("time") + scale_linetype_discrete(name = "Quantity")
# print(g)
start <- 10
g <- ggplot(subset(state1, Time >= start), 
            aes(x = Time, y = quantile95, colour = Method, linetype = Run))
g <- g + geom_line(alpha = 1)
g <- g + geom_line(aes(x = Time, y = TRUESTATES[start:T,1]), colour = "darkgreen", alpha = 0.2,
                   linetype = 1)
# pdf(file = "FirstStatePrediction.pdf", height = 5, width = 10)
print(g)
# dev.off()


sum(rep(TRUESTATES[,1], nrunsbsmc) > subset(state1, Method == "BSMC")[,"quantile95"]) / nrunsbsmc
sum(rep(TRUESTATES[,1], nrunssmc2) > subset(state1, Method == "SMC2")[,"quantile95"]) / nrunssmc2
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 1))$quantile95)
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 2))$quantile95)
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 3))$quantile95)
# sum(TRUESTATES[,1] > (subset(state1, Method == "SMC2" & Run == 4))$quantile95)



