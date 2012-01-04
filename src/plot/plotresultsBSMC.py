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


class PlotResultsBSMC(PlotResults):
    def __init__(self, resultsfolder, RDatafile):
        self.method = "BSMC"
        self.color = "orange"
        PlotResults.__init__(self, resultsfolder, RDatafile)
        self.Rcode += """pdf(file = pdffile, useDingbats = FALSE, title = "%s results")\n""" % self.method
        self.parametersHaveBeenLoaded = False
    def everything(self):
        self.addEvidence()
        self.allParameters()
        self.addObservations()
        self.addPredictedObs()
        self.close()
    def singleparameter(self, parameterindex):
        self.histogramparameter(parameterindex)
    def loadparameters(self):
        self.Rcode += \
"""
indexhistory <- length(savingtimes)
t <- savingtimes[indexhistory]
thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
"""
        self.parametersHaveBeenLoaded = True
    def histogramparameter(self, parameterindex):
        if not(self.parametersHaveBeenLoaded):
            self.loadparameters()
        self.Rcode += \
"""
i <- %(parameterindex)i
g <- qplot(x = thetas[,i], geom = "blank") + geom_histogram(aes(y = ..density..)) + geom_density(fill = "%(color)s", alpha = 0.5)
g <- g + xlab(%(parametername)s)
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

