# The file requires two RData files to be available. They should be
# runs from SMC^2 and BSMC on the same model, same observations, same horizon T.
rm(list = ls())
gc()
graphics.off()
resultsfile <- "/home/pierre/Dropbox/py-smc2/results/BSMC-SVonefactor-synthetic-T1000-N250000-h0.100(0).RData"
load(file = resultsfile)
BSMCpred1 <- predictedstate1
BSMCpred2 <- predictedstate2
BSMCpredSqObs <- predictedsquaredobs
keepfromBSMC <- c("keepfromBSMC", "BSMCpred1", "BSMCpred2", "BSMCpredSqObs", "observations")
rm(list = ls()[!(ls() %in% keepfromBSMC)])
resultsfile <- "/home/pierre/Dropbox/py-smc2/results/SMC2-SVonefactor-synthetic-T1000-independent-dynamicNx250-Ntheta1000(0).RData"
load(file = resultsfile)
SMC2pred1 <- predictedstate1
SMC2pred2 <- predictedstate2
SMC2predSqObs <- predictedsquaredobs
rm(list = ls()[!(ls() %in% c(keepfromBSMC, "SMC2pred1", "SMC2pred1", "SMC2predSqObs", "truestates",
                             "T"))])

library(ggplot2)
library(gridExtra)

g <- qplot(x = 1:T, y = truestates[,1], geom = "line", colour = "True states")
g <- g + geom_line(aes(y = BSMCpred1[,1], colour = "BSMC mean"))
g <- g + geom_line(aes(y = BSMCpred1[,2], colour = "conf interval"))
g <- g + geom_line(aes(y = BSMCpred1[,3], colour = "conf interval"))
g <- g + ylim(0, 5) + xlim(250, 1000) + scale_colour_discrete(name = "")
g <- g + xlab("time") + ylab(expression(v[t]))
# print(g)

g2 <- qplot(x = 1:T, y = truestates[,1], geom = "line", colour = "True states")
g2 <- g2 + geom_line(aes(y = SMC2pred1[,1], colour = "SMC2 mean"))
g2 <- g2 + geom_line(aes(y = SMC2pred1[,2], colour = "conf interval"))
g2 <- g2 + geom_line(aes(y = SMC2pred1[,3], colour = "conf interval"))
g2 <- g2 + ylim(0, 5) + xlim(250, 1000) + scale_colour_discrete(name = "")
g2 <- g2 + xlab("time") + ylab(expression(v[t]))
# print(g2)
sqrt(sum((SMC2pred1[,1] - truestates[,1])**2))
sqrt(sum((BSMCpred1[,1] - truestates[,1])**2))

pdf(file = "/home/pierre/compareSVpred.pdf")
grid.arrange(g, g2)
start <- 250
g <- qplot(x = start:T, y = observations[start:T]**2, geom = "line", colour = "squared observations")
g <- g + geom_line(aes(y = SMC2predSqObs[start:T,1], colour = "SMC2 predicted mean"))
g <- g + geom_line(aes(y = SMC2predSqObs[start:T,2], colour = "SMC2 90 pct confidence"))
g <- g + geom_line(aes(y = SMC2predSqObs[start:T,3], colour = "SMC2 90 pct confidence"))
g <- g + scale_colour_discrete(name = "")
g <- g + xlab("time") + ylab("squared observations")
g2 <- qplot(x = start:T, y = observations[start:T]**2, geom = "line", colour = "squared observations")
g2 <- g2 + geom_line(aes(y = BSMCpredSqObs[start:T,1], colour = "BSMC predicted mean"))
g2 <- g2 + geom_line(aes(y = BSMCpredSqObs[start:T,2], colour = "BSMC 90 pct confidence"))
g2 <- g2 + geom_line(aes(y = BSMCpredSqObs[start:T,3], colour = "BSMC 90 pct confidence"))
g2 <- g2 + scale_colour_discrete(name = "")
g2 <- g2 + xlab("time") + ylab("squared observations")
grid.arrange(g, g2)
dev.off()


