
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/periodic/synthetic/"
library(ggplot2)
setwd(resultsfolder)
theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15))

### Load SMC2 indep results
SMC2Indresultsfile <- paste(resultsfolder, "SMC2-T250-independent-dynamicNx250-Ntheta1000(", sep = "")
SMC2IndDF <- data.frame()
Nruns <- 2
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
### Load Liu and West results
BSMCresultsfile <- paste(resultsfolder, "BSMC-T250-N500000(", sep = "")
#Nruns <- 4
BSMCDF <- data.frame()
for (run in 1:Nruns){
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
#BSMCDF <- SMC2IndDF

smc2Indweightedmeans <- matrix(nrow = Nruns, ncol = nbparameters)
BSMCDFweightedmeans <- matrix(nrow = Nruns, ncol = nbparameters)

for (run in 1:Nruns){
    subrunSMC2Ind <- subset(SMC2IndDF, as.double(Run) == run)
    subrunBSMC <- subset(BSMCDF, as.double(Run) == run)
    for (indexparam in 1:nbparameters){
        smc2Indweightedmeans[run, indexparam] <- weighted.mean(x = 
          subrunSMC2Ind[[indexparam]], weights = subrunSMC2Ind$w)
        BSMCDFweightedmeans[run, indexparam] <- weighted.mean(x = 
          subrunBSMC[[indexparam]], weights = subrunBSMC$w)
    }
}


compar <- rbind(data.frame(smc2Indweightedmeans, method = "smc2 independent"), 
                data.frame(BSMCDFweightedmeans, method = "liu and west"))
names(compar) <- c(paste("Theta", 1:nbparameters, sep = ""), "method")
compar <- melt(compar, id = c("method"))
compar$method <- as.factor(compar$method)
g <- ggplot(compar, aes(x = variable, y = value, colour = method))
g <- g + geom_boxplot()
#g <- g + geom_boxplot() + scale_y_log()
print(g)

variances <- rbind(data.frame(var = rbind(apply(smc2Indweightedmeans, 2, var)), method = "smc2 independent"),
data.frame(var = rbind(apply(BSMCDFweightedmeans, 2, var)), method = "liu and west"))
names(variances) <- c(paste("Theta", 1:nbparameters, sep = ""), "method")
variances <- melt(variances, id = c("method"))
variances$method <- as.factor(variances$method)
g <- ggplot(variances, aes(x = variable, y = value, fill= method))
g <- g + geom_bar(position = "dodge") + scale_y_log()
X11();print(g)

