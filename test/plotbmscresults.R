rm(list = ls())
setwd("~/Dropbox/py-smc2/")
load("results.RData")

dim(thetahistory)
par(mfrow = c(2, 1))
hist(thetahistory[1, 1,])
hist(thetahistory[1, 2,])
indexhistory <- length(savingtimes)
t <- savingtimes[indexhistory]
# w <- weighthistory[indexhistory,]
# w <- w / sum(w)
thetasDF <- as.data.frame(t(thetahistory[indexhistory,,]))
# thetasDF <- cbind(thetas, w)
names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""))
head(thetasDF)