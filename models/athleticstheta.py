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

from numpy import random, power, sqrt, exp, zeros_like, zeros, \
        ones, mean, average, prod, log
from src.parametermodel import ParameterModel

def logdprior(parameters, hyperparameters):
    """ takes transformed parameters """
    nu_part = parameters[0] - hyperparameters["nu_rate"] * exp(parameters[0])
    xi_part = parameters[1] - hyperparameters["xi_rate"] * exp(parameters[1])
    sigma_part = parameters[2] - hyperparameters["sigma_rate"] * exp(parameters[2])
    return nu_part + xi_part + sigma_part

def rprior(size, hyperparameters):
    """ returns untransformed parameters """
    nu = random.exponential(scale = 1 / hyperparameters["nu_rate"], size = size)
    xi = random.exponential(scale = 1 / hyperparameters["xi_rate"], size = size)
    sigma = random.exponential(scale = 1 / hyperparameters["sigma_rate"], size = size)
    parameters = zeros((3, size))
    parameters[0, :] = nu
    parameters[1, :] = xi
    parameters[2, :] = sigma
    return parameters

## Prior hyperparameters
hyperparameters = { \
        "nu_rate": 0.2, \
        "xi_rate": 0.5, \
        "sigma_rate": 0.2}

modeltheta = ParameterModel(name = "Athletics records", dimension = 3)
modeltheta.setHyperparameters(hyperparameters)
modeltheta.setPriorlogdensity(logdprior)
modeltheta.setPriorgenerator(rprior)
modeltheta.setParameterNames(["expression(nu)", "expression(xi)", "expression(sigma)"])
modeltheta.setTransformation(["log", "log", "log"])

Rpriornu = "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["nu_rate"]
Rpriorxi = "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["xi_rate"]
Rpriorsigma = "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["sigma_rate"]
modeltheta.setRprior([Rpriornu, Rpriorxi, Rpriorsigma])



#modeltheta.plottingInstructions = ["No observations", "No evidence"]
#modeltheta.additionalPlots = """
#filteringDF <- data.frame(observations)
#filteringDF$year <- 1976:2010
#filteringDF$xfiltered <- filteredfirststate
#filteringDF$CIlow <- filteredCIlow
#filteringDF$CIhigh <- filteredCIhigh
#names(filteringDF) <- c("y1", "y2", "year", "x", "CIlow", "CIhigh")
#g <- ggplot(filteringDF, aes(x = year))
#g <- g + geom_line(aes(y = x), size = 1, linetype = 2, colour = "blue")
#g <- g + geom_line(aes(y = CIlow), size = 1, colour = "blue")
#g <- g + geom_line(aes(y = CIhigh), size = 1, colour = "blue")
#g <- g + geom_point(aes(y = y1), size = 4, colour = "green")
#g <- g + geom_point(aes(y = y2), size = 4, colour = "red")
#g <- g + ylab("Time (seconds)") + ylim(450, 540)
#print(g)
#smoothingDF <- data.frame(observations)
#smoothingDF$year <- 1976:2010
#smoothingDF$xsmoothed <- smoothedmeansfirststate35
#smoothingDF$CIlow <- smoothedmeansCIlow35
#smoothingDF$CIhigh <- smoothedmeansCIhigh35
#names(smoothingDF) <- c("y1", "y2", "year", "x", "CIlow", "CIhigh")
#g <- ggplot(smoothingDF, aes(x = year))
#g <- g + geom_line(aes(y = x), size = 1, linetype = 2, colour = "black")
#g <- g + geom_line(aes(y = CIlow), size = 1, colour = "blue")
#g <- g + geom_line(aes(y = CIhigh), size = 1, colour = "blue")
#g <- g + geom_point(aes(y = y1), size = 4, colour = "green")
#g <- g + geom_point(aes(y = y2), size = 4, colour = "red")
#g <- g + ylab("Time (seconds)")
#g <- g + ylim(450, 540)
#print(g)
#"""
#
#InverseGammaTemplate = """
#priorfunction <- function(x){
#    shape <- %.5f 
#    scale <- %.5f
#    return(scale**shape / gamma(shape) * x**(- shape - 1) * exp(-scale / x))
#}
#"""
#Rpriornu =  InverseGammaTemplate % (hyperparameters["nu_shape"], hyperparameters["nu_scale"])
#Rpriorxi = """priorfunction <- function(x) dnorm(x, mean = %.5f, sd = %.5f)""" \
#    % (hyperparameters["xi_mean"], hyperparameters["xi_sd"])

