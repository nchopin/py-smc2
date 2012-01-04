rm(list = ls())
gc()
graphics.off()
library(ggplot2)
setwd("/home/pierre/Dropbox/py-smc2/finalplots/")
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/finalplots/SVSP500-results/"
singlefactorfile <- paste(resultsfolder, "SVonefactor/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
twofactorfile <- paste(resultsfolder, "SVfixedrho/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
fullfile <- paste(resultsfolder, "SVmultifactor/SP500recent/SMC2-T753-independent-dynamicNx100-Ntheta2000(", sep = "")
nrunssmc2 <- 2
#load(file = paste(SMC2Indresultsfile, 0, ").RData", sep = ""))
#NBPARAM <- dim(thetahistory)[2]
#OBSERVATIONS <- observations

evid1fact<- loadSomeQuantity(basename = singlefactorfile, nruns = nrunssmc2,
                                  somequantity = "evidences", somename = "LogEvid",
                                  functiontoapply = function(x) cumsum(log(x)))
evid1fact$Model <- "SV one factor"
evid2fact <- loadSomeQuantity(basename = twofactorfile, nruns = nrunssmc2,
                                  somequantity = "evidences", somename = "LogEvid",
                                  functiontoapply = function(x) cumsum(log(x)))
evid2fact$Model <- "SV multifactor"
evidFull <- loadSomeQuantity(basename = fullfile, nruns = nrunssmc2,
                                  somequantity = "evidences", somename = "LogEvid",
                                  functiontoapply = function(x) cumsum(log(x)))
evidFull$Model <- "SV full"

allevid <- rbind(evid1fact, evid2fact, evidFull)
allevid$Run <- factor(allevid$Run)
allevid$Model<- factor(allevid$Model)

evid <- subset(allevid, Run == 1)
evid$Run <- factor(evid$Run)
evid2versus1 <- subset(evid, Model == "SV multifactor")$LogEvid - subset(evid, Model == "SV one factor")$LogEvid
evidFullversus1 <- subset(evid, Model == "SV full")$LogEvid - subset(evid, Model == "SV one factor")$LogEvid

otherevid <- subset(allevid, Run == 2)
otherevid$Run <- factor(otherevid$Run)
otherevid2versus1 <- subset(otherevid, Model == "SV multifactor")$LogEvid - 
                     subset(otherevid, Model == "SV one factor")$LogEvid
otherevidFullversus1 <- subset(otherevid, Model == "SV full")$LogEvid - 
                        subset(otherevid, Model == "SV one factor")$LogEvid

g <- qplot(x = 1:753, y = evid2versus1, geom = "line", colour = "Multifactor without leverage")
g <- g + geom_line(aes(y = evidFullversus1, colour = "Multifactor with leverage"))
g <- g + scale_colour_discrete(name = "Model:", h = c(0, 0), c = c(0, 0), l = c(0, 50))
g <- g + xlab("Time") + ylab("Log evidence comparison")
g <- g + opts(legend.position = c(0.4, 0.8), legend.background = theme_rect(fill = "white"),
              legend.title = theme_text(size = 25, hjust = -0.2), legend.text = theme_text(size = 20))
pdf(file = "SP500evidencecomparison.pdf")
print(g)
dev.off()

#g <- qplot(x = 1:753, y = otherevid2versus1, geom = "line", colour = "Multifactor without leverage")
#g <- g + geom_line(aes(y = otherevidFullversus1, colour = "Multifactor with leverage"))
#g <- g + scale_colour_discrete(name = "Model:")
#g <- g + xlab("Time") + ylab("Log evidence compared to the single factor model")
#g <- g + opts(legend.position = c(0.4, 0.8), legend.background = theme_rect(fill = "white"),
#              legend.title = theme_text(size = 20, hjust = -0.2), legend.text = theme_text(size = 15))
#X11(); print(g)



