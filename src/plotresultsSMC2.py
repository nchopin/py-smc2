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
from src.plotresults import PlotResults

class PlotResultsSMC2(PlotResults):
    def __init__(self, resultsfolder, RDatafile):
        self.method = "SMC2"
        self.color = "blue"
        PlotResults.__init__(self, resultsfolder, RDatafile)
        self.Rcode += """pdf(file = pdffile, useDingbats = FALSE, title = "%s results")\n""" % self.method
        self.parametersHaveBeenLoaded = False

    def everything(self):
        self.acceptancerate()
        self.ESS()
        self.addComputingTime()
        self.evidence()
        self.allParameters()
        self.addObservations()
        self.close()

    def singleparameter(self, parameterindex):
        self.histogramparameter(parameterindex)

    def loadparameters(self):
        self.Rcode += \
"""
indexhistory <- length(savingtimes)
t <- savingtimes[indexhistory]
w <- weighthistory[indexhistory,]
w <- w / sum(w)
#thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
#thetasDF <- cbind(thetas, w)
#names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w")
"""
        self.parametersHaveBeenLoaded = True
    def acceptancerate(self):
        self.Rcode += \
"""
acceptratiodataframe <- as.data.frame(cbind(resamplingindices, acceptratios))
g <- ggplot(data = acceptratiodataframe, aes(x = resamplingindices, y= acceptratios))
g <- g + geom_point(size = 4) + geom_line() + xlab("iterations") + ylab("acceptance rates")
g <- g + xlim(0, T) + ylim(0, 1) 
print(g)
"""
    def ESS(self):
        self.Rcode += \
"""
Ntheta <- dim(thetahistory)[3]
ESSdataframe <- as.data.frame(cbind(1:length(ESS), ESS))
g <- ggplot(data = ESSdataframe, aes(x = V1, y= ESS))
g <- g + geom_line() + xlab("iterations") + ylab("ESS") + ylim(0, Ntheta)
print(g)
"""
    def evidence(self):
        if not("No evidence" in self.plottingInstructions):
            self.Rcode += \
"""
evidencedataframe <- as.data.frame(cbind(1:length(evidences), evidences))
g <- ggplot(data = evidencedataframe, aes(x = V1, y= evidences))
g <- g + geom_line() + xlab("iterations") + ylab("evidence")
print(g)
"""

    def histogramparameter(self, parameterindex):
        if not(self.parametersHaveBeenLoaded):
            self.loadparameters()
        self.Rcode += \
"""
indexhistory <- length(savingtimes)
w <- weighthistory[indexhistory,]
w <- w / sum(w)
i <- %(parameterindex)i
g <- qplot(x = thetahistory[indexhistory,i,], weight = w, geom = "blank") + 
  geom_histogram(aes(y = ..density..)) + geom_density(fill = "blue", alpha = 0.5) +
    xlab(%(parametername)s)
""" % {"parameterindex": parameterindex + 1, "parametername": self.parameternames[parameterindex], "color": self.color}
        if hasattr(self.modeltheta, "truevalues"):
            self.Rcode += \
"""
g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
"""
        if hasattr(self.modeltheta, "Rfunctionlist"):
            self.Rcode += \
"""
%s
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
""" % self.modeltheta.Rfunctionlist[parameterindex]
            if hasattr(self.modelx, "Rlikelihood"):
                self.Rcode += \
"""
%s
trueposterior <- function(x) priorfunction(x) * truelikelihood(x)
g <- g + stat_function(fun = trueposterior, colour = "green", size = 2)
""" % self.modelx.Rlikelihood[parameterindex]
        self.Rcode += \
"""
print(g)
"""
    def addComputingTime(self):
        self.Rcode += \
"""
g <- qplot(x = 1:T, y = cumsum(computingtimes), geom = "line",
           ylab = "computing time", xlab = "iteration")
print(g)
"""

