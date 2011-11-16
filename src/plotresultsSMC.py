###################################################
#    This file is part of py-smc2.
#    http://code.google.com/p/py-smc2/
#
#    py-smc2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    py-smc2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with py-smc2.  If not, see <http://www.gnu.org/licenses/>.
###################################################

#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import os, os.path
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, cov, isnan, zeros_like, \
        var, isinf, linalg, pi, dot

class PlotResultsSMC:
    def __init__(self, RDatafile, resultsfolder = ""):
        self.RDatafile = RDatafile
        if len(resultsfolder) == 0:
            self.resultsfolder = os.path.dirname(self.RDatafile) 
        else:
            self.resultsfolder = resultsfolder
        self.pdffile = os.path.basename(self.RDatafile.replace(".RData", ".pdf"))
        self.plotresultsfile = self.pdffile.replace(".pdf", "-plots.R")
        self.plotresultsfile = os.path.join(self.resultsfolder, self.plotresultsfile)
        self.plottingInstructions = []
        self.Rcode = \
"""
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "%(resultspath)s/"
resultsfile <- "%(RDatafile)s"
pdffile <- "%(pdffile)s"
library(ggplot2)
setwd(resultsfolder)
cat("working with file", resultsfile, "...\n")
load(file = resultsfile)
pdf(file = pdffile, useDingbats = FALSE, title = "SMC filtering results")
""" % {"RDatafile": self.RDatafile, "resultspath": self.resultsfolder, \
        "pdffile": self.pdffile}
        self.addFilteringMean()
        self.addTrajectories()
        #self.addComputingTime()
        #self.addObservations()
        #self.addTrueStates()
    def addObservations(self):
        if not("No observations" in self.plottingInstructions):
            self.Rcode += \
"""
observationsDF <- cbind(data.frame(observations), 1:length(observations))
names(observationsDF) <- c("y", "index")
g <- ggplot(data = observationsDF, aes(x = index, y = y)) 
g <- g + geom_line() + ylab("observations")
print(g)
"""
    def addTrueStates(self):
        self.Rcode += \
"""
truestatesDF <- cbind(data.frame(truestates), 1:length(truestates))
names(truestatesDF) <- c(paste("x", 1:statedimension, sep =""), "index")
g <- ggplot(data = melt(truestatesDF, id = "index"), aes(x = index, y = value, colour = variable))
g <- g + geom_line() + ylab("true states")
print(g)
"""
    def addFilteringMean(self):
        self.Rcode += \
"""
#statesDF <- data.frame(x = truestates, xhat = meanpath, index = 1:length(truestates))
#statesDF <- melt(statesDF, id = "index")
#g <- ggplot(data = statesDF, aes(x = index, y = value, colour = variable)) 
#g <- g + geom_line() + ylab("states")
g <- qplot(x = 1:length(truestates), y = meanpath, colour = "SMC mean", geom = "line") +
    geom_line(aes(y = truestates, colour = "Truth")) +
    xlab("time") + ylab("hidden state") + scale_colour_discrete(name = "")
print(g)
"""
    def addTrajectories(self):
        self.Rcode += \
"""
if (exists("trajectories")){
    library(foreach)
    df <- foreach(indexdim = 1:statedimension, .combine = rbind) %do% {
      trajectoriesDF  <- data.frame(x = truestates[,indexdim], index = 1:length(truestates[,1]))
      for (index in 1:dim(trajectories)[3]){
        trajectoriesDF <- cbind(trajectoriesDF, trajectories[2:(T+1),indexdim,index])
      }
      names(trajectoriesDF) <- c("x", "index", paste("trajectory", 1:dim(trajectories)[3], sep = ""))
      trajectoriesDF <- melt(trajectoriesDF, id = "index")
      trajectoriesDF$dim <- paste("X", indexdim, sep = "")
      trajectoriesDF
    }
    g <- ggplot(data = df, aes(x = index, y = value, colour = variable))
    g <- g + facet_wrap( ~ dim)
    g <- g + geom_line() + ylab("states") + scale_colour_manual(values = c("blue", colors()[grep("orange", colors())]))
    g <- g + opts(legend.position = "none")
    print(g)
}
"""
    def addComputingTime(self):
        self.Rcode += \
"""
g <- qplot(x = 1:T, y = cumsum(computingtimes), geom = "line",
           ylab = "computing time", xlab = "iteration")
print(g)
"""
    def setDLM(self, string):
        self.Rcode += string
        self.Rcode += \
"""
KF <- function(observations, somedlm){
  # the notation comes from the package dlm
  T <- length(observations)
  m <- rep(0, T + 1); C <- rep(1, T + 1)
  a <- rep(0, T); R <- rep(0, T)
  f <- rep(0, T); Q <- rep(0, T)
  m[1] <- somedlm$m0; C[1] <- somedlm$C0
  for (t in 1:T){
    a[t] <- somedlm$GG * m[t]
    R[t] <- somedlm$GG * C[t] * somedlm$GG + somedlm$W
    f[t] <- somedlm$FF * a[t]
    Q[t] <- somedlm$FF * R[t] * somedlm$FF + somedlm$V
    m[t+1] <- a[t] + R[t] * somedlm$FF * (1 / Q[t]) * (observations[t] - f[t])
    C[t+1] <- R[t] - R[t] * somedlm$FF * (1 / Q[t]) * somedlm$FF * R[t]
  }
  return(list(observations = observations, NextObsMean = f, NextObsVar = Q,
              NextStateMean = a, NextStatevar = R,
              FiltStateMean = m[2:(T+1)], FiltStateVar = C[2:(T+1)]))
}
"""
    def addKalmanComparison(self):
        self.Rcode += \
"""
kalmanresults <- KF(observations, dlm)
#xhat <- apply(X=xhistory[,1,2:(T+1)], MARGIN = 2, FUN=mean)
g <- qplot(x = 1:T, y = meanpath, geom = "line", colour = "SMC mean") +
     geom_line(aes(y = kalmanresults$FiltStateMean, colour = "KF mean"), alpha = 0.) + 
     geom_point(aes(y = kalmanresults$FiltStateMean, colour = "KF mean")) +
     geom_point(aes(y = truestates, colour= "Truth")) + 
     xlab("time") + ylab("hidden state") + scale_colour_discrete(name = "")
print(g)
if (exists("xhistory")){
    for (index in 1:length(savingtimes)){
        sometime <- savingtimes[index]
        if (sometime > 0){
            truefilteringdistribution <- function(x, sometime){
              dnorm(x, mean = kalmanresults$FiltStateMean[sometime], 
                        sd = sqrt(kalmanresults$FiltStateVar[sometime]))
            }
            truedist <- function(x) truefilteringdistribution(x, sometime)
            g <- qplot(x = xhistory[,,index], geom = "blank") + 
              aes_string(x = paste("xhistory[,,", index, "]")) +
              stat_function(fun = truedist, colour = "red") + 
              geom_density(fill = "blue", alpha = 0.3) + xlab(paste("hidden state at time", sometime)) + ylab("density")
            print(g)
        }
    }
}
"""
    def close(self):
        self.Rcode += \
"""
dev.off()
"""
        f = open(self.plotresultsfile, "w")
        f.write(self.Rcode)
        f.close()





