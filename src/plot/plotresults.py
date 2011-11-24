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

class PlotResults:
    def __init__(self, resultsfolder, RDatafile):
        self.resultsfolder = resultsfolder
        self.RDatafile = RDatafile
        self.pdffile = self.RDatafile.replace(".RData", ".pdf")
        self.plotresultsfile = self.pdffile.replace(".pdf", "-plots.R")
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
""" % {"RDatafile": self.RDatafile, "resultspath": self.resultsfolder, \
        "pdffile": self.pdffile}
    def setModelTheta(self, modeltheta):
        self.modeltheta = modeltheta
        #self.additionalPlots = self.modeltheta.additionalPlots
        self.Rcode += \
"""
nbparameters <- %i
""" % self.modeltheta.parameterdimension
    def setModelX(self, modelx):
        self.modelx = modelx
    def setParameterNames(self, names):
        self.parameterdimension = len(names)
        self.parameternames = names
    def allParameters(self):
        for parameterindex in range(self.parameterdimension):
            self.singleparameter(parameterindex)
    def addObservations(self):
        self.Rcode += \
"""
observationsDF <- cbind(data.frame(observations), 1:length(observations))
names(observationsDF) <- c("y", "index")
g <- ggplot(data = observationsDF, aes(x = index, y = y)) 
g <- g + geom_line() + ylab("observations")
print(g)
"""
    def addPredictedObs(self):
            self.Rcode += \
"""
if (exists("truestates") && is.null(dim(truestates))){
    truestates <- as.matrix(truestates, ncol = 1)
}
predictedquantities <- grep(patter="predicted", x = ls(), value = TRUE)
if (T > 25){
    start <- 10
} else {
    start <- 1
}
for (name in predictedquantities){
    ystr <- paste(name, "[start:T,1]", sep = "")
    yqt1 <- paste(name, "[start:T,2]", sep = "")
    yqt2 <- paste(name, "[start:T,3]", sep = "")
    g <- qplot(x = start:T, geom = "blank")
    g <- g + geom_line(aes_string(y = ystr, colour = paste("'", name, "'", sep = "")))
    g <- g + geom_line(aes_string(y = yqt1, colour = paste("'", name, "quantile'", sep = "")))
    g <- g + geom_line(aes_string(y = yqt2, colour = paste("'", name, "quantile'", sep = "")))
    if (name == "predictedstate1" && exists("truestates")){
            g <- g + geom_line(aes(y = truestates[start:T,1], colour = "True states"))
    } else {
        if (name == "predictedstate2" && exists("truestates")){
            g <- g + geom_line(aes(y = truestates[start:T,2], colour = "True states"))
        } else {
            if (name == "predictedobservations"){
                g <- g + geom_line(aes(y = observations[start:T,1], colour = "observations"))
            }
            if (name == "predictedsquaredobs"){
                g <- g + geom_line(aes(y = observations[start:T,1]**2, colour = "squared observations"))
            }
        }
    }
    g <- g + opts(legend.position = "bottom", legend.direction = "horizontal")
    g <- g + xlab("time") + ylab(name) + scale_colour_discrete(name = "")
    print(g)
}
"""
    def addEvidence(self):
        self.Rcode += \
"""
evidencedataframe <- as.data.frame(cbind(1:length(evidences), evidences))
g <- ggplot(data = evidencedataframe, aes(x = V1, y= evidences))
g <- g + geom_line() + xlab("iterations") + ylab("evidence")
print(g)
"""

    def close(self):
        #self.Rcode += self.additionalPlots
        self.Rcode += \
"""
dev.off()
"""
        f = open(self.plotresultsfile, "w")
        f.write(self.Rcode)
        f.close()




