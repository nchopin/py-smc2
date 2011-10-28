
rm(list = ls())
gc()
graphics.off()
source("/home/pierre/Dropbox/py-smc2/rgraphs/plotutils.R")

resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/SVonefactor/synthetic/"
setwd(resultsfolder)

# load true parameters and nbparameters from one of the results file
adPMCMCresultsfile <- paste(resultsfolder, "adPMCMC-T", 1000, "-Iter50000-Nx500(", sep = "")
load(file = paste(adPMCMCresultsfile, 0, ").RData", sep = ""))

### load results
## SMC^2
TIMES <- c(250, 500, 1000)

nrunssmc2 <- 3
SMC2Indresultsfile <- paste(resultsfolder, "SMC2-T1000-independent-dynamicNx250-Ntheta1000(", sep = "")
#smc2 <- loadSMCresult(basename = SMC2Indresultsfile, nruns = nrunssmc2, time = TIME)

## BSMC
nrunsbsmc <- 3
BSMCsmooth <- "-h0.010"
BSMCresultsfile <- paste(resultsfolder, "BSMC-T1000-N250000", BSMCsmooth, "(", sep = "")
#bsmc <- loadSMCresult(basename = BSMCresultsfile, nruns = nrunsbsmc, time = TIME)

library(foreach)
smc2alltimes <- foreach (timeindex = 1:length(TIMES), .combine = rbind) %do%{
    smc2 <- loadSMCresult(basename = SMC2Indresultsfile, nruns = nrunssmc2, time = TIMES[timeindex])
    smc2
}
bsmcalltimes <- foreach (timeindex = 1:length(TIMES), .combine = rbind) %do%{
    bsmc <- loadSMCresult(basename = BSMCresultsfile, nruns = nrunsbsmc, time = TIMES[timeindex])
    bsmc 
}

smc2alltimes$method = "smc2"
bsmcalltimes$method = "liu and west"

## Adaptive PMCMC
nrunsadpmcmc <- 1
adpmcmcalltimes <- foreach (timeindex = 1:length(TIMES), .combine = rbind) %do%{
    adPMCMCresultsfile <- paste(resultsfolder, "adPMCMC-T", TIMES[timeindex], "-Iter50000-Nx500(", sep = "")
    adpmcmc <- loadMCMCresult(basename = adPMCMCresultsfile, nruns = nrunsadpmcmc, 
                              time = TIMES[timeindex], burnin = 10000)
    adpmcmc
}
adpmcmcalltimes$method = "PMCMC"


all <- rbind(smc2alltimes, bsmcalltimes, adpmcmcalltimes)


i <- 4
g <- ggplot(all, 
            aes_string(x = paste("Theta", i, sep = ""), weight = "w", colour = "method", fill = "Run"))
g <- g + facet_wrap(method ~ Time)
g <- g + geom_density(alpha = 0.1)
g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
g <- g + xlim(0, 0.5) + ylab("Density")
#print(g)

ggsave(g, filename = paste("Theta", i, ".pdf", sep = ""))

### Adaptive PMCMC
#nrunsadpmcmc <- 1
#adPMCMCresultsfile <- paste(resultsfolder, "adPMCMC-T", TIME, "-Iter50000-Nx500(", sep = "")
#adpmcmc <- loadMCMCresult(basename = adPMCMCresultsfile, nruns = nrunsadpmcmc, time = TIME, burnin = 10000)
#
## load true parameters and nbparameters from one of the results file
#load(file = paste(adPMCMCresultsfile, 0, ").RData", sep = ""))
#
#adpmcmc$method = "adpmcmc"
#smc2$method = "smc2"
#bsmc$method = "liu and west"
#all <- rbind(adpmcmc, smc2, bsmc)

### 
i <- 4
g <- ggplot(all, aes_string(x = paste("Theta", i, sep = ""), weight = "w", linetype = "Run", colour = "method"))
g <- g + geom_density(alpha = 0.1)
g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
g <- g + xlim(0, 0.5) + ylab("Density")
print(g)
ggsave(plot = g, filename = paste("theta", i, "-T", TIME, ".pdf", sep = ""))

i <- 5
g <- ggplot(all, aes_string(x = paste("Theta", i, sep = ""), weight = "w", fill = "Run", colour = "method"))
g <- g + geom_density(alpha = 0.1)
g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
g <- g + xlim(0, 0.03) + ylab("Density")
print(g)
ggsave(plot = g, filename = paste("theta", i, "-T", TIME, ".pdf", sep = ""))


adpmcmcmeans <- matrix(nrow = nrunsadpmcmc, ncol = nbparameters)
smc2means <- matrix(nrow = nrunssmc2, ncol = nbparameters)
bsmcmeans <- matrix(nrow = nrunsbsmc, ncol = nbparameters)
for (run in 1:nrunsadpmcmc){
    subadPMCMC <- subset(adpmcmc, as.double(Run) == run)
    for (indexparam in 1:nbparameters){
        adpmcmcmeans[run, indexparam] <- mean(subadPMCMC[[indexparam]])
    }
}
for (run in 1:nrunsbsmc){
    subrunBSMC <- subset(bsmc, as.double(Run) == run)
    for (indexparam in 1:nbparameters){
        bsmcmeans[run, indexparam] <- weighted.mean(x = 
          subrunBSMC[[indexparam]], weights = subrunBSMC$w)
    }
}
for (run in 1:nrunssmc2){
    subrunSMC2Ind <- subset(smc2, as.double(Run) == run)
    for (indexparam in 1:nbparameters){
        smc2means[run, indexparam] <- weighted.mean(x = 
          subrunSMC2Ind[[indexparam]], weights = subrunSMC2Ind$w)
    }
}
compar <- rbind(data.frame(smc2means, method = "smc2 independent"), 
                data.frame(bsmcmeans, method = "liu and west"),
                data.frame(adpmcmcmeans, method = "adaptive PMCMC"))
names(compar) <- c(paste("Theta", 1:nbparameters, sep = ""), "method")
compar <- melt(compar, id = c("method"))
compar$method <- as.factor(compar$method)
g <- ggplot(compar, aes(x = variable, y = value, colour = method))
g <- g + geom_boxplot()
X11();print(g)

comparLast <- subset(compar, variable %in% c("Theta4"))
comparLast$variable <- factor(comparLast$variable)
g <- ggplot(comparLast, aes(x = variable, y = value, colour = method))
g <- g + geom_point()
X11();print(g)

#variances <- rbind(data.frame(var = rbind(apply(smc2Indweightedmeans, 2, var)), method = "smc2 independent"),
#data.frame(var = rbind(apply(BSMCDFweightedmeans, 2, var)), method = "liu and west"))
#names(variances) <- c(paste("Theta", 1:nbparameters, sep = ""), "method")
#variances <- melt(variances, id = c("method"))
#variances$method <- as.factor(variances$method)
#g <- ggplot(variances, aes(x = variable, y = value, fill= method))
#g <- g + geom_bar(position = "dodge") + scale_y_log()
#X11();print(g)


#i <- 4
#g <- ggplot(SMC2IndDF, aes(SMC2IndDF[[i]], weight = w, fill = Run))
#g <- g + geom_density(alpha = 0.2)
#g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
#g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
#g <- g + opts(legend.position = "none") + ylab("Density")
#print(g)
#g <- ggplot(BSMCDF, aes(BSMCDF[[i]], weight = w, fill = Run))
#g <- g + geom_density(alpha = 0.2)
#g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
#g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
#g <- g + opts(legend.position = "none") + ylab("Density")
#print(g)
#

# head(SMC2IndDF)
# pdf(file = paste(resultsfolder, "compareparameters.pdf"))
# for (i in 1:nbparameters){
#     g <- ggplot(SMC2IndDF, aes_string(x = paste("Theta", i, sep = ""), weight = "w", fill = "Run"))
#     g <- g + geom_density(alpha = 0.2)
#     g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
#     g <- g + opts(legend.position = "none") + ylab("Density")
#     print(g)
# }
# dev.off()

#adPMCMCDF2 <- adPMCMCDF[, -8]

