
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/pySMC2/results/SVonefactor"
resultsfile <- "/home/pierre/Dropbox/pySMC2/results/SVonefactor/synthetic/SMC2-T1000-dynamicNx50-Nth1000("
library(ggplot2)
setwd(resultsfolder)
theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15),
strip.text.x = theme_text(size=20))
trueparameters <- c(0, 0, 0.5, 0.0625, 0.01)

allrunsDF <- data.frame()
for (savingtime in 1:4){
    actualtime <- paste("T = ", savingtime * 250, sep = "")
    for (run in 1:10){
        load(file = paste(resultsfile, run - 1, ").RData", sep = ""))
        indexhistory <- savingtime
        w <- weighthistory[indexhistory,]
        w <- w / sum(w)
        thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
        thetasDF <- cbind(thetas, w, run, actualtime)
        names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run", "Time")
        allrunsDF <- rbind(allrunsDF, thetasDF)
    }
}
allrunsDF$Run <- as.factor(allrunsDF$Run)
allrunsDF$Time <- as.factor(allrunsDF$Time)

expressionlabels <- c(expression(mu), expression(beta), expression(xi), expression(omega^2), expression(lambda))
    i <- 1
    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
    g <- g + geom_density(alpha = 0.2) + xlab(expressionlabels[i])
    priorfunction <- function(x) dnorm(x, sd = 1.41421)
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
    g <- g + opts(legend.position = "none") + ylab("Density")
    g <- g + xlim(0, 0.05)
    g <- g + facet_grid(. ~ Time)
    X11(width = 20, height = 5); print(g)
    ggsave(file = "SVonefactor-mu-Concentration.pdf", useDingbats = FALSE)
    i <- 2
    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
    g <- g + geom_density(alpha = 0.2) + xlab(expressionlabels[i])
    priorfunction <- function(x) dnorm(x, sd = 1.41421)
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
    g <- g + opts(legend.position = "none") + ylab("Density")
    g <- g + xlim(0, 0.05)
    g <- g + facet_grid(. ~ Time)
    X11(width = 20, height = 5); print(g)
    ggsave(file = "SVonefactor-beta-Concentration.pdf", useDingbats = FALSE)
    i <- 3
    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
    g <- g + geom_density(alpha = 0.2) + xlab(expressionlabels[i])
    priorfunction <- function(x) dnorm(x, sd = 1.41421)
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
    g <- g + opts(legend.position = "none") + ylab("Density")
    g <- g + xlim(0, 0.05)
    g <- g + facet_grid(. ~ Time)
    X11(width = 20, height = 5); print(g)
    ggsave(file = "SVonefactor-xi-Concentration.pdf", useDingbats = FALSE)
    i <- 4
    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
    g <- g + geom_density(alpha = 0.2) + xlab(expressionlabels[i])
    priorfunction <- function(x) dnorm(x, sd = 1.41421)
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
    g <- g + opts(legend.position = "none") + ylab("Density")
    g <- g + xlim(0, 0.05)
    g <- g + facet_grid(. ~ Time)
    X11(width = 20, height = 5); print(g)
    ggsave(file = "SVonefactor-omega2-Concentration.pdf", useDingbats = FALSE)
    i <- 5
    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
    g <- g + geom_density(alpha = 0.2) + xlab(expressionlabels[i])
    priorfunction <- function(x) dnorm(x, sd = 1.41421)
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
    g <- g + opts(legend.position = "none") + ylab("Density")
    g <- g + xlim(0, 0.05)
    g <- g + facet_grid(. ~ Time)
    X11(width = 20, height = 5); print(g)
    ggsave(file = "SVonefactor-lambda-Concentration.pdf", useDingbats = FALSE)






    i <- 5
    #TIME <- 250
    #g <- ggplot(subset(allrunsDF, Run == 1), aes(subset(allrunsDF, Run == 1)[[i]], weight = w, fill = Time))
    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
    g <- g + geom_density(alpha = 0.2) + xlab(expressionlabels[i])
#g <- g + geom_histogram(aes(y=..density..), colour = "black", fill = "white")
    priorfunction <- function(x) dnorm(x, sd = 1.41421)
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
    g <- g + opts(legend.position = "none") + ylab("Density")
    g <- g + xlim(0, 0.1)
    g <- g + facet_grid(. ~ Time)
    X11(width = 10, height = 5); print(g)
    ggsave(file = "SVonefactor-lambda-Concentration.pdf", useDingbats = FALSE)

    i <- 2
    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
    g <- g + geom_density(alpha = 0.2) + xlab(expression(beta))
    priorfunction <- function(x) dnorm(x, sd = 1.41421)
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
    g <- g + opts(legend.position = "none") + ylab("Density")
    print(g)
    ggsave(file = paste("SVonefactor-beta-T", actualtime, ".pdf", sep = ""), useDingbats = FALSE)

#    i <- 3
#    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
#    g <- g + geom_density(alpha = 0.2) + xlab(expression(xi))
#    priorfunction <- function(x) dexp(x, rate = 0.20000)
#    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
#    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
#    g <- g + opts(legend.position = "none") + ylab("Density")
#    print(g)
#    ggsave(file = "SVonefactor-xi.pdf", useDingbats = FALSE)
#
#    i <- 4
#    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
#    g <- g + geom_density(alpha = 0.2) + xlab(expression(omega^2))
#    priorfunction <- function(x) dexp(x, rate = 0.20000)
#    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
##g <- g + scale_x_log10()
#    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1, n = 2000)
#    g <- g + opts(legend.position = "none") + ylab("Density")
#    g <- g + xlim(0, 2)
#    print(g)
#    ggsave(file = "SVonefactor-omega2.pdf", useDingbats = FALSE)
#
#    i <- 5
#    g <- ggplot(allrunsDF, aes(allrunsDF[[i]], weight = w, fill = Run))
#    g <- g + geom_density(alpha = 0.2) + xlab(expression(lambda))
#    priorfunction <- function(x) dexp(x, rate = 1.00000)
#    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
#    g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
##g <- g + scale_x_log10()
#    g <- g + xlim(0, 0.05)
#    g <- g + opts(legend.position = "none") + ylab("Density")
#    print(g)
#    ggsave(file = "SVonefactor-lambda.pdf", useDingbats = FALSE)



