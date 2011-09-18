
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/"
library(ggplot2)
setwd(resultsfolder)
theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15))

### Load SMC2 indep results
SMC2Indresultsfile <- "/home/pierre/Dropbox/py-smc2/results/ind("
SMC2IndDF <- data.frame()
Nruns <- 10
for (run in 1:Nruns){
    load(file = paste(SMC2Indresultsfile, run - 1, ").RData", sep = ""))
    indexhistory <- length(savingtimes)
    t <- savingtimes[indexhistory]
    w <- weighthistory[indexhistory,]
    w <- w / sum(w)
    thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
    thetasDF <- cbind(thetas, w, run)
    names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run")
    SMC2IndDF <- rbind(SMC2IndDF, thetasDF)
}
SMC2IndDF$Run <- as.factor(SMC2IndDF$Run)
### Load SMC2 RW results
SMC2RWresultsfile <- "/home/pierre/Dropbox/py-smc2/results/thetalogistic/synthetic/SMC2-T1000-dynamicNx250-Nth1000("
Nruns <- 10
SMC2RWDF <- data.frame()
for (run in 1:Nruns){
    load(file = paste(SMC2RWresultsfile, run - 1, ").RData", sep = ""))
    indexhistory <- length(savingtimes)
    t <- savingtimes[indexhistory]
    w <- weighthistory[indexhistory,]
    w <- w / sum(w)
    thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
    thetasDF <- cbind(thetas, w, run)
    names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run")
    SMC2RWDF <- rbind(SMC2RWDF, thetasDF)
}
SMC2RWDF$Run <- as.factor(SMC2RWDF$Run)
### Load Liu and West results
BSMCresultsfile <- "/home/pierre/Dropbox/py-smc2/results/thetalogistic/synthetic/BSMC-T1000-N500000("
Nruns <- 9
BSMCDF <- data.frame()
for (run in 3:Nruns){
    load(file = paste(BSMCresultsfile, run - 1, ").RData", sep = ""))
    indexhistory <- length(savingtimes)
    t <- savingtimes[indexhistory]
    w <- weighthistory[indexhistory,]
    w <- w / sum(w)
    thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
    thetasDF <- cbind(thetas, w, run)
    names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run")
    BSMCDF <- rbind(BSMCDF, thetasDF)
}
BSMCDF$Run <- as.factor(BSMCDF$Run)



smc2Indweightedmeans <- matrix(nrow = Nruns, ncol = nbparameters)
smc2RWweightedmeans <- matrix(nrow = Nruns, ncol = nbparameters)
smc2RWweightedmeans <- matrix(nrow = Nruns, ncol = nbparameters)

for (run in 1:Nruns){
    subrunSMC2IndDF <- subset(SMC2IndDF, as.double(Run) == run)
    subrunSMC2RWDF <- subset(SMC2RWDF, as.double(Run) == run)
    for (indexparam in 1:nbparameters){
        smc2Indweightedmeans[run, indexparam] <- weighted.mean(x = 
          subrunSMC2IndDF[[indexparam]], weights = subrunSMC2IndDF$w)
        smc2RWweightedmeans[run, indexparam] <- weighted.mean(x = 
          subrunSMC2RWDF[[indexparam]], weights = subrunSMC2RWDF$w)
    }
}

apply(smc2Indweightedmeans, 2, var)
apply(smc2RWweightedmeans, 2, var)

compar <- rbind(data.frame(smc2Indweightedmeans, method = "smc2 independent"), 
                data.frame(smc2RWweightedmeans, method = "smc2 random walk"))
names(compar) <- c(paste("Theta", 1:nbparameters, sep = ""), "method")
compar <- melt(compar, id = c("method"))
compar$method <- as.factor(compar$method)
g <- ggplot(compar, aes(x = variable, y = value, colour = method))
g <- g + geom_boxplot()
print(g)


