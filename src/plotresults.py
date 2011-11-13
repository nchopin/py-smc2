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
""" % {"RDatafile": self.RDatafile, "resultspath": self.resultsfolder, \
        "pdffile": self.pdffile}
    def setModelTheta(self, modeltheta):
        self.modeltheta = modeltheta
        if hasattr(self.modeltheta, "truevalues"):
            self.Rcode += \
"""
truevalues <- c(%s)
""" % ",".join([str(value) for value in self.modeltheta.truevalues])
        if len(self.modeltheta.plottingInstructions) > 0:
            self.plottingInstructions = self.modeltheta.plottingInstructions
        self.additionalPlots = self.modeltheta.additionalPlots
    def setModelX(self, modelx):
        self.modelx = modelx
        if hasattr(self.modelx, "RLinearGaussian"):
            self.Rcode += "\n" + self.modelx.RLinearGaussian
    def setParameters(self, names):
        self.parameterdimension = len(names)
        self.parameternames = names
    def allParameters(self):
        for parameterindex in range(self.parameterdimension):
            self.singleparameter(parameterindex)
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

    def close(self):
        self.Rcode += self.additionalPlots
        self.Rcode += \
"""
dev.off()
"""
        f = open(self.plotresultsfile, "w")
        f.write(self.Rcode)
        f.close()




