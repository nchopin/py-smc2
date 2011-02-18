
rm(list = ls())
gc()
library(ggplot2)
theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15))
#theme_set(newtheme)
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/pySMC2/results/"
resultsfile <- "/home/pierre/Dropbox/pySMC2/results/SMC2-athletics-athletics-best-two-T35-dynamicNx2000-Nth5000("
setwd(resultsfolder)

#pdf(file = "results3runs.pdf", useDingbats = FALSE, title = "SMC2 results")
allrunsDF <- data.frame()
allNx <- data.frame()
allAccept <- data.frame()
allbeatprobaDF <- data.frame()
allsmoothvalues <- data.frame()
for (run in 1:10){
    load(file = paste(resultsfile, run - 1, ").RData", sep = ""))
    Nxlist2 <- rep(Nxlist, 1, each = 2)
    increaseindices2 <- c(0, rep(increaseindices[2:length(increaseindices)], 1, each = 2), T)
    allNx <- rbind(allNx, cbind(Nxlist2, increaseindices2, run))
    allAccept <- rbind(allAccept, cbind(resamplingindices, acceptratios, run))
    indexhistory <- length(savingtimes)
    t <- savingtimes[indexhistory]
    w <- weighthistory[indexhistory,]
    w <- w / sum(w)
    thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
    thetasDF <- cbind(thetas, w, run, smoothedvaluesbeat199335, smoothedvaluesbeat198535)
    names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run", "smooth1993", "smooth1985")
    allrunsDF <- rbind(allrunsDF, thetasDF)
    beatprobaDF <- data.frame(observations)
    beatprobaDF$year <- 1976:2010
    beatprobaDF$smoothbeat1985 <- smoothedmeansbeat198535
    beatprobaDF$filterbeat1985 <- filteredbeat1985
    beatprobaDF$smoothbeat1993 <- smoothedmeansbeat199335
    beatprobaDF$filterbeat1993 <- filteredbeat1993
    names(beatprobaDF) <- c("y1", "y2", "year", "smooth1985", "filter1985", "smooth1993", "filter1993")
    beatprobaDF$Run <- run
    allbeatprobaDF <- rbind(allbeatprobaDF, beatprobaDF)
}

allrunsDF$Run <- as.factor(allrunsDF$Run)
names(allNx) <- c("Nx", "iterations", "Run")
allNx$Run <- as.factor(allNx$Run)
names(allAccept) <- c("iterations", "acceptratios", "Run")
allAccept$Run <- as.factor(allAccept$Run)
allbeatprobaDF$Run <- as.factor(allbeatprobaDF$Run)


g <- ggplot(data = allAccept, aes(x = iterations, y= acceptratios, colour = Run))
g <- g + geom_point(size = 4) + geom_line() + xlab("iterations") + ylab("acceptance rates")
g <- g + xlim(0, T) + ylim(0, 1) 
print(g)

i <- 1
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(nu))
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
g <- g + opts(axis.title.x=theme_text(size=20), axis.title.y=theme_text(size=20, angle = 90))
g <- g + opts(axis.text.x=theme_text(size=15), axis.text.y=theme_text(size=15))
g <- g + xlim(0, 5)
pdf(file = "athletics-nu.pdf", useDingbats = FALSE)
print(g)
dev.off()
i <- 2
g <- ggplot(allrunsDF, aes(-allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(xi))
priorfunction <- function(x) dexp(-x, rate = 0.50000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
g <- g + xlim(-0.8, 0)
pdf(file = "athletics-xi.pdf", useDingbats = FALSE)
print(g)
dev.off()
i <- 3
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(sigma))
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
g <- g + xlim(2, 7)
pdf(file = "athletics-sigma.pdf", useDingbats = FALSE)
print(g)
dev.off()

names(allbeatprobaDF)
g <- ggplot(allbeatprobaDF, aes(group = year, x = year))
g <- g + geom_boxplot(aes(y = smooth1985), colour = "black") + stat_summary(fun.y = mean, geom = "line", aes(group = 1, y = smooth1985), colour = "black")
g <- g + geom_boxplot(aes(y = smooth1993), colour = "red") + stat_summary(fun.y = mean, geom = "line", aes(group = 1, y = smooth1993), colour = "red")
g <- g + geom_vline(xintercept = 1993, linetype = 3, size = 1)
g <- g + scale_y_log10() + ylab("Probability") + xlab("Year")
pdf(file = "athletics-probaofbeating.pdf", useDingbats = FALSE)
print(g)
dev.off()

g <- ggplot(allbeatprobaDF, aes(group = year, x = year))
g <- g + geom_boxplot(aes(y = smooth1985), colour = "black") + stat_summary(fun.y = mean, geom = "line", aes(group = 1, y = smooth1985), colour = "black")
g <- g + geom_boxplot(aes(y = smooth1993), colour = "red") + stat_summary(fun.y = mean, geom = "line", aes(group = 1, y = smooth1993), colour = "red")
g <- g + geom_boxplot(aes(y = smooth1993 / smooth1985)) + stat_summary(fun.y = mean, geom = "line", aes(group = 1, y = smooth1993 / smooth1985), linetype = 2)
g <- g + geom_vline(xintercept = 1993, linetype = 3, size = 1)
g <- g + scale_y_log10() + ylab("Probability") + xlab("Year")
pdf(file = "athletics-probaofbeatingWithConditional.pdf", useDingbats = FALSE)
print(g)
dev.off()

g <- ggplot(allbeatprobaDF, aes(group = year, x = year))
g <- g + geom_boxplot(aes(y = smooth1993 / smooth1985)) + stat_summary(fun.y = mean, geom = "line", aes(group = 1, y = smooth1993 / smooth1985), linetype = 2)
g <- g + geom_vline(xintercept = 1993, linetype = 3, size = 1)
g <- g + scale_y_log10() + ylab("Probability") + xlab("Year")
pdf(file = "athletics-probaConditional.pdf", useDingbats = FALSE)
print(g)
dev.off()

atYear1993 <- subset(allbeatprobaDF, year == 1993)
probaOfInterest <- atYear1993$smooth1993 / atYear1993$smooth1985
print("probability of interest:")
print(probaOfInterest)
print(summary(probaOfInterest))
print(sd(probaOfInterest))

#g <- ggplot(allrunsDF, aes(x = allrunsDF[[2]], y = smooth1985, weight = w))
#g <- g + geom_point(aes(alpha = w, size = w), colour = "blue", shape = 1)
#g <- g + geom_point(aes(alpha = w, size = w, y = smooth1993), colour = "red", shape = 2)
#g <- g + opts(legend.position = "none")
#g <- g + xlab(expression(sigma)) + ylab("probability")
#g <- g + scale_y_log()
#print(g)
#
#g <- ggplot(allbeatprobaDF, aes(group = year, x = year))
#g <- g + geom_boxplot(aes(y = smooth1993 / smooth1985)) + stat_summary(fun.y = mean, geom = "line", aes(group = 1, y = smooth1993 / smooth1985), linetype = 2)
#g <- g + ylab("probability")
#pdf(file = "athletics-probaofinterest.pdf", useDingbats = FALSE)
#print(g)
#dev.off()


