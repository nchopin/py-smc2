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
from src.plot.plotresults import PlotResults


class PlotResultsPMMH(PlotResults):
    def __init__(self, resultsfolder, RDatafile, burnin):
        self.method = "adPMCMC"
        self.color = "red"
        PlotResults.__init__(self, resultsfolder, RDatafile)
        self.Rcode += """pdf(file = pdffile, useDingbats = FALSE, title = "%s results")\n""" % self.method
        self.parametersHaveBeenLoaded = False
        self.burnin = burnin

    def everything(self):
        self.allParameters()
        self.addObservations()
        self.close()

    def singleparameter(self, parameterindex):
        self.histogramparameter(parameterindex)
        self.traceplot(parameterindex)
    def loadparameters(self):
        self.Rcode += \
"""
thetasDF <- data.frame(t(chain))
names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""))
thetasDF$iterations <- 1:(dim(thetasDF)[1])
thetasDF <- subset(thetasDF, iterations > %(burnin)i)
""" % {"burnin": self.burnin}
        self.parametersHaveBeenLoaded = True

    def histogramparameter(self, parameterindex):
        if not(self.parametersHaveBeenLoaded):
            self.loadparameters()
        self.Rcode += \
"""
i <- %(parameterindex)i
g <- ggplot(thetasDF, aes(thetasDF[[i]]))
g <- g + geom_histogram(aes(y=..density..), colour = "black", fill = "white")
g <- g + geom_density(fill = "%(color)s", alpha = 0.5) + xlab(%(parametername)s)
""" % {"parameterindex": parameterindex + 1, "parametername": self.parameternames[parameterindex], "color": self.color}
        self.Rcode += \
"""
if (exists("trueparameters")){
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
}
"""
        if hasattr(self.modeltheta, "Rprior"):
            self.Rcode += \
"""
%s
g <- g + stat_function(fun = priorfunction, colour = "red", linetype = 1, size = 1)
""" % self.modeltheta.Rprior[parameterindex]
        self.Rcode += \
"""
print(g)
"""

    def traceplot(self, parameterindex):
        if not(self.parametersHaveBeenLoaded):
            self.loadparameters()
        self.Rcode += \
"""
i <- %(parameterindex)i + 1
g <- ggplot(thetasDF, aes(x = iterations, y = thetasDF[[i]]))
g <- g + geom_line() + ylab(%(parametername)s)
print(g)
""" % {"parameterindex": parameterindex, "parametername": self.parameternames[parameterindex], "color": self.color}

