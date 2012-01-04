rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

explabels <- c(expression(nu), expression(xi), expression(sigma))

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/Athletics-results/"
SMC2Indresultsfile <- paste(resultsfolder, "SMC2-T35-independent-dynamicNx250-Ntheta1000(", sep = "")
nrunssmc2 <- 10
load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
NBPARAM <- dim(thetahistory)[2]
OBSERVATIONS <- observations

library(foreach)
df <- foreach(run = 1:nrunssmc2, .combine = rbind) %do% {
    load(file = paste(SMC2Indresultsfile, run - 1, ").RData", sep = ""))
    quantityDF <- as.data.frame(cbind(smoothedmeansbeat198535, 
                                      smoothedmeansbeat199335, 1976:2010, run))
    names(quantityDF) <- c("beat1985", "beat1993", "Time", "Run")
    quantityDF
}
df$Run <- as.factor(df$Run)
df$Cond <- df$beat1993 / df$beat1985
head(df)
cat("mean cond prob at time 1993", mean(subset(df, Time == "1993")$Cond), "\n")
cat("sd cond prob at time 1993", sd(subset(df, Time == "1993")$Cond), "\n")

g <- ggplot(df, aes(group = Time, x = Time))
g <- g + geom_boxplot(aes(y = beat1985), colour = "black")
g <- g + stat_summary(fun.y = mean, geom = "line", 
                      aes(group = 1, y = beat1985), colour = "black")
g <- g + geom_boxplot(aes(y = beat1993), colour = "red")
g <- g + stat_summary(fun.y = mean, geom = "line", 
                      aes(group = 1, y = beat1993), colour = "red")
g <- g + geom_boxplot(aes(y = Cond), colour = "blue")
g <- g + stat_summary(fun.y = mean, geom = "line", 
                      aes(group = 1, y = Cond), colour = "blue")
g <- g + scale_y_log10() + ylab("Probability") + xlab("Year")
g <- g + geom_vline(xintercept = 1993, linetype = 3)
pdf(file = "Athl-boxplots.pdf")
print(g)
dev.off()