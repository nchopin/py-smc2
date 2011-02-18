
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/pySMC2/results/"
resultsfileSOPF <- "/home/pierre/Dropbox/pySMC2/results/SOPF-SVonefactor-synthetic-T100-N10000("
resultsfileSMC2 <- "/home/pierre/Dropbox/pySMC2/results/SMC2-SVonefactor-synthetic-T100-dynamicNx50-Nth100("
library(ggplot2)
setwd(resultsfolder)
trueparameters <- c(0, 0, 0.5, 0.0625, 0.01)

sopfDF <- data.frame()
for (run in 1:50){
    load(file = paste(resultsfileSOPF, run - 1, ").RData", sep = ""))
    indexhistory <- length(savingtimes)
    finaltime <- savingtimes[indexhistory]
    particles <- allreducedparticles[[indexhistory]]
    nbparticles = dim(particles)[1]
    cat("Run", run, ", number of particles:", nbparticles, "\n")
    w <- allcounts[[indexhistory]]
    thetas <- as.data.frame(particles)
    # thetasDF <- cbind(thetas, w)
    thetasDF <- cbind(thetas, w, run)
    names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run")
    sopfDF <- rbind(sopfDF, thetasDF)
}

smc2DF <- data.frame()
allNx <- data.frame()
allAccept <- data.frame()
for (run in 1:49){
    load(file = paste(resultsfileSMC2, run - 1, ").RData", sep = ""))
    Nxlist2 <- rep(Nxlist, 1, each = 2)
    increaseindices2 <- c(0, rep(increaseindices[2:length(increaseindices)], 1, each = 2), T)
    allNx <- rbind(allNx, cbind(Nxlist2, increaseindices2, run))
    allAccept <- rbind(allAccept, cbind(resamplingindices, acceptratios, run))
    indexhistory <- length(savingtimes)
    t <- savingtimes[indexhistory]
    w <- weighthistory[indexhistory,]
    w <- w / sum(w)
    thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
    thetasDF <- cbind(thetas, w, run)
    names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run")
    smc2DF <- rbind(smc2DF, thetasDF)
}

sopfweightedmeans <- matrix(nrow = 10, ncol = 5)
smc2weightedmeans <- matrix(nrow = 10, ncol = 5)
for (run in 1:50){
    subrunsopfDF <- subset(sopfDF, as.double(Run) == run)
    subrunsmc2DF <- subset(smc2DF, as.double(Run) == run)
    for (indexparam in 1:5){
        sopfweightedmeans[run, indexparam] <- weighted.mean(x = subrunsopfDF[[indexparam]], weights = subrunsopfDF$w)
        smc2weightedmeans[run, indexparam] <- weighted.mean(x = subrunsmc2DF[[indexparam]], weights = subrunsmc2DF$w)
    }
}

par(mfrow = c(1, 2))
boxplot(sopfweightedmeans, ylim = c(-5, 10))
boxplot(smc2weightedmeans, ylim = c(-5, 10))
compar <- data.frame(sopfweightedmeans)
compar$method = "sopf"

compar2 <- data.frame(smc2weightedmeans)
compar2$method = "smc2"

compar <- data.frame(rbind(cbind(sopfweightedmeans, "sopf"), cbind(smc2weightedmeans, "smc2")))
names(compar) <- c(paste("Theta", 1:nbparameters, sep = ""), "method")
compar <- melt(compar, id = c("method"))
compar$method <- as.factor(compar$method)
compar$value <- as.double(compar$value)
g <- ggplot(compar, aes(x = variable, y = value, colour = method))
g <- g + geom_boxplot()
print(g)


i <- 1
g <- ggplot(thetasDF, aes(thetasDF[[i]], weight = w))  
g <- g + geom_histogram(aes(y=..density..), colour = "black", fill = "white")
g <- g + geom_density(fill = "green", alpha = 0.5) + xlab(expression(mu))
print(g)

i <- 2
g <- ggplot(thetasDF, aes(thetasDF[[i]], weight = w))  
g <- g + geom_histogram(aes(y=..density..), colour = "black", fill = "white")
g <- g + geom_density(fill = "green", alpha = 0.5) + xlab(expression(beta))
print(g)

i <- 3
g <- ggplot(thetasDF, aes(thetasDF[[i]], weight = w))  
g <- g + geom_histogram(aes(y=..density..), colour = "black", fill = "white")
g <- g + geom_density(fill = "green", alpha = 0.5) + xlab(expression(xi))
print(g)

i <- 4
g <- ggplot(thetasDF, aes(thetasDF[[i]], weight = w))  
g <- g + geom_histogram(aes(y=..density..), colour = "black", fill = "white")
g <- g + geom_density(fill = "green", alpha = 0.5) + xlab(expression(omega^2))
print(g)

i <- 5
g <- ggplot(thetasDF, aes(thetasDF[[i]], weight = w))  
g <- g + geom_histogram(aes(y=..density..), colour = "black", fill = "white")
g <- g + geom_density(fill = "green", alpha = 0.5) + xlab(expression(lambda))
print(g)

observationsDF <- cbind(data.frame(observations), 1:length(observations))
names(observationsDF) <- c("y", "index")
g <- ggplot(data = observationsDF, aes(x = index, y = y)) 
g <- g + geom_line() + ylab("observations")
print(g)

observationsDF <- cbind(data.frame(observations**2), 1:length(observations))
names(observationsDF) <- c("y2", "index")
g <- ggplot(data = observationsDF, aes(x = index, y = y2)) 
g <- g + geom_line() + ylab("squared observations")
print(g)

dev.off()
