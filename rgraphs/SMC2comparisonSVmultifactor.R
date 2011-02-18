
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/pySMC2/results/SVfixedrho/"
resultsfile <- "/home/pierre/Dropbox/pySMC2/results/SVfixedrho/SP500easy/SMC2-T753-dynamicNx125-Nth1000("
library(ggplot2)
setwd(resultsfolder)
theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15))

#trueparameters <- c(0, 0, 0.5, 0.0625, 0.01)

allrunsDF <- data.frame()
allNx <- data.frame()
allAccept <- data.frame()
for (run in 1:2){
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
    thetasDF <- cbind(thetas, w, run)
    names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run")
    allrunsDF <- rbind(allrunsDF, thetasDF)
}

allrunsDF$Run <- as.factor(allrunsDF$Run)
names(allNx) <- c("Nx", "iterations", "Run")
allNx$Run <- as.factor(allNx$Run)
names(allAccept) <- c("iterations", "acceptratios", "Run")
allAccept$Run <- as.factor(allAccept$Run)

pdf(file = "SVfixedrho-SP500easy.pdf", pointsize = 24, useDingbats = FALSE)

g <- ggplot(data = allAccept, aes(x = iterations, y= acceptratios, colour = Run))
g <- g + geom_point(size = 4) + geom_line() + xlab("Iterations") + ylab("Acceptance rates")
g <- g + xlim(0, T) + ylim(0, 1) 
g <- g + opts(legend.position = "none")
print(g)
#ggsave(file = "SVmultifactor-acceptrates.pdf", useDingbats = FALSE)

g <- ggplot(data = allNx, aes(x = iterations, y = Nx, colour = Run))
g <- g + geom_line() + geom_point()
g <- g + opts(legend.position = "none") + xlab("Iterations")
print(g)
#ggsave(file = "SVmultifactor-Nx.pdf", useDingbats = FALSE)

i <- 1
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(mu))
priorfunction <- function(x) dnorm(x, sd = 1.41421)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
print(g)
#ggsave(file = "SVmultifactor-mu.pdf", useDingbats = FALSE)


i <- 2
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(beta))
priorfunction <- function(x) dnorm(x, sd = 1.41421)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
print(g)
#ggsave(file = "SVmultifactor-beta.pdf", useDingbats = FALSE)

i <- 3
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(xi))
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
print(g)
#ggsave(file = "SVmultifactor-xi.pdf", useDingbats = FALSE)

i <- 4
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(omega^2))
priorfunction <- function(x) dexp(x, rate = 0.20000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1, n = 2000)
g <- g + opts(legend.position = "none") + ylab("Density")
print(g)
#ggsave(file = "SVmultifactor-omega2.pdf", useDingbats = FALSE)
7
i <- 5
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(lambda[1]))
priorfunction <- function(x) dexp(x, rate = 1.00000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
print(g)

i <- 6
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(lambda[2] - lambda[1]))
priorfunction <- function(x) dexp(x, rate = 0.50000)
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
print(g)

i <- 7
g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
g <- g + geom_density(alpha = 0.2) + xlab(expression(w[1]))
priorfunction <- function(x) 1
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
g <- g + opts(legend.position = "none") + ylab("Density")
print(g)

observationsDF <- cbind(data.frame(observations), 1:length(observations))
names(observationsDF) <- c("y", "index")
g <- ggplot(data = observationsDF, aes(x = index, y = y)) 
g <- g + geom_line() + ylab("Observations") + xlab("Time")
print(g)
#ggsave(file = "SVmultifactor-observations.pdf", useDingbats = FALSE)

observationsDF <- cbind(data.frame(observations**2), 1:length(observations))
names(observationsDF) <- c("y2", "index")
g <- ggplot(data = observationsDF, aes(x = index, y = y2)) 
g <- g + geom_line() + ylab("Squared observations") + xlab("Time")
#g <- g + geom_vline(xintercept = allAccept$iterations[allAccept$Run == 1], linetype = 2, size = 0.5, colour = "red")
print(g)
#ggsave(file = "SVmultifactor-observations2.pdf", useDingbats = FALSE)

dev.off()

