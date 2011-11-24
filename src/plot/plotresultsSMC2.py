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
        self.addEvidence()
        self.allParameters()
        self.addObservations()
        self.addPredictedObs()
        self.addFiltered()
        self.close()

    def singleparameter(self, parameterindex):
        self.histogramparameter(parameterindex)

    def loadparameters(self):
        self.Rcode += \
"""
#indexhistory <- length(savingtimes)
#t <- savingtimes[indexhistory]
#w <- weighthistory[indexhistory,]
#w <- w / sum(w)
"""
        self.parametersHaveBeenLoaded = True
    def acceptancerate(self):
        self.Rcode += \
"""
g <- qplot(x = resamplingindices, y = acceptratios, geom = "line", colour = "acceptance rates")
#g <- g + geom_line(aes(y = guessAR, colour = "guessed acceptance rates"))
g <- g + xlim(0, T) + ylim(0, 1) +  opts(legend.position = "bottom", legend.direction = "horizontal")
g <- g + scale_colour_discrete(name = "") + xlab("time") + ylab("acceptance rates")
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
#g <- ggplot(data = ESSdataframe, aes(x = V1, y= ESS))
#g <- g + geom_line() + xlab("iterations (square root scale)") + ylab("ESS (log)") + ylim(0, Ntheta)
#g <- g + scale_x_sqrt() + scale_y_log()
#print(g)
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
g <- qplot(x = thetahistory[indexhistory,i,], weight = w, geom = "blank")
g <- g + geom_histogram(aes(y = ..density..)) + geom_density(fill = "blue", alpha = 0.5)
g <- g + xlab(%(parametername)s)
if (exists("trueparameters")){
    g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
}
g <- g + opts(legend.position = "bottom", legend.direction = "horizontal")
""" % {"parameterindex": parameterindex + 1, "parametername": self.parameternames[parameterindex], "color": self.color}
        if hasattr(self.modeltheta, "Rprior"):
            self.Rcode += \
"""
%s
g <- g + stat_function(fun = priorfunction, aes(colour = "prior"), linetype = 1, size = 1)
g <- g + scale_colour_discrete(name = "")
""" % self.modeltheta.Rprior[parameterindex]
        self.Rcode += \
"""
print(g)
"""
    def addComputingTime(self):
        self.Rcode += \
"""
g <- qplot(x = 1:T, y = cumsum(computingtimes), geom = "line",
           ylab = "computing time (square root scale)", xlab = "iteration")
g <- g + scale_y_sqrt()
print(g)
"""
    def addFiltered(self):
        self.Rcode += \
"""
if (exists("truestates") && is.null(dim(truestates))){
    truestates <- as.matrix(truestates, ncol = 1)
}
filteredquantities <- grep(patter="filtered", x = ls(), value = TRUE)
for (name in filteredquantities){
    g <- qplot(x = 1:T, geom = "blank") 
    g <- g + geom_line(aes_string(y = name,
    colour = paste("'",name, "'", sep = "")))
    if (name == "filteredstate1" && exists("truestates")){
            g <- g + geom_line(aes(y = truestates[,1], colour = "True states"))
    } else {
        if (name == "filteredstate2" && exists("truestates"))
            g <- g + geom_line(aes(y = truestates[,2], colour = "True states"))
    }
    g <- g + xlab("time") + ylab(name) + scale_colour_discrete(name = "")
    g <- g + opts(legend.position = "bottom", legend.direction = "horizontal")
    print(g)
}
"""













