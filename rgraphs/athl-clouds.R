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


## posterior distribution
smc2thetas <- loadSMCThetaParticles(basename = SMC2Indresultsfile, 
                                   nruns = nrunssmc2, time = 35)


library(foreach)
df <- foreach(run = 1:nrunssmc2, .combine = rbind) %do% {
    load(file = paste(SMC2Indresultsfile, run - 1, ").RData", sep = ""))
    quantityDF <- as.data.frame(cbind(smoothedvaluesbeat198535, 
                                      smoothedvaluesbeat199335, smoothedvaluesweights, 
                                      run))
    names(quantityDF) <- c("beat1985", "beat1993", "updatedlogw", "Run")
    quantityDF
}

df$updatedw <- exp(df$updatedlogw)
df$Run <- factor(df$Run)
for (r in levels(df$Run)){
    S <- sum(df[df$Run == r,"updatedw"])
    df[df$Run == r,"updatedw"] <-df[df$Run == r,"updatedw"] / S 
}
allrunsDF <- df
allrunsDF$Theta1 <- smc2thetas$Theta1
allrunsDF$Theta2 <- smc2thetas$Theta2
allrunsDF$Theta3 <- smc2thetas$Theta3

head(allrunsDF)

g <- ggplot(allrunsDF, aes(x = Theta1))
g <- g + geom_point(aes(y = beat1985, alpha = updatedw, 
                        shape = "1985", colour = "1985"))
g <- g + geom_point(aes(y = beat1993, alpha = updatedw, 
                        shape = "1993", colour = "1993"))
g <- g + scale_colour_discrete(name = "", 
                                      h = c(0, 200), c = c(100, 100), l = c(50, 20))
g <- g + opts(legend.position = "none")
g <- g + xlab(expression(nu)) + ylab("Probability") + scale_y_log()
# print(g)
ggsave(g, filename="Athl-CloudTheta1.png", width = 7, height = 7, dpi = 150)

g <- ggplot(allrunsDF, aes(x = -Theta2))
g <- g + geom_point(aes(y = beat1985, alpha = updatedw, 
                        shape = "1985", colour = "1985"))
g <- g + geom_point(aes(y = beat1993, alpha = updatedw, 
                        shape = "1993", colour = "1993"))
g <- g + scale_colour_discrete(name = "", 
                                      h = c(0, 200), c = c(100, 100), l = c(50, 20))
g <- g + opts(legend.position = "none")
g <- g + xlab(expression(xi)) + ylab("Probability") + scale_y_log()
ggsave(g, filename="Athl-CloudTheta2.png", width = 7, height = 7, dpi = 150)

g <- ggplot(allrunsDF, aes(x = Theta3))
g <- g + geom_point(aes(y = beat1985, alpha = updatedw, 
                        shape = "1985", colour = "1985"))
g <- g + geom_point(aes(y = beat1993, alpha = updatedw, 
                        shape = "1993", colour = "1993"))
g <- g + scale_colour_discrete(name = "", 
                                      h = c(0, 200), c = c(100, 100), l = c(50, 20))
g <- g + opts(legend.position = "none")
g <- g + xlab(expression(sigma)) + ylab("Probability") + scale_y_log()
ggsave(g, filename="Athl-CloudTheta3.png", width = 7, height = 7, dpi = 150)




