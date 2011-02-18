#! /usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import random, power, sqrt, exp, zeros_like, zeros, \
        ones, mean, average, prod, log
from scipy.stats import norm, truncnorm, gamma
from src.models import ParameterModel


#### See src/models.py for explanations about the model functions.
#### most functions here take transformed parameters, except rprior

#def invgamma_logpdf(logx, shape, scale):
#    return (-shape - 1) * logx - scale / exp(logx)

def logdprior(parameters, hyperparameters):
    """ takes transformed parameters """
    #nu_part = parameters[0] + invgamma_logpdf(parameters[0], hyperparameters["nu_shape"], hyperparameters["nu_scale"])
    nu_part = parameters[0] - hyperparameters["nu_rate"] * exp(parameters[0])
    #xi_part = - 0.5 / (hyperparameters["xi_sd"]**2) * (parameters[1] - hyperparameters["xi_mean"])**2
    xi_part = parameters[1] - hyperparameters["xi_rate"] * exp(parameters[1])
    #sigma_part = parameters[2] + invgamma_logpdf(parameters[2], hyperparameters["sigma_shape"], hyperparameters["sigma_scale"])
    sigma_part = parameters[2] - hyperparameters["sigma_rate"] * exp(parameters[2])
    return nu_part + xi_part + sigma_part

def rprior(size, hyperparameters):
    """ returns untransformed parameters """
    #nu = 1 / gamma.rvs(hyperparameters["nu_shape"], scale = hyperparameters["nu_scale"], size = size)
    #xi = norm.rvs(size = size, loc = hyperparameters["xi_mean"], scale = hyperparameters["xi_sd"])
    nu = random.exponential(scale = 1 / hyperparameters["nu_rate"], size = size)
    xi = random.exponential(scale = 1 / hyperparameters["xi_rate"], size = size)
    sigma = random.exponential(scale = 1 / hyperparameters["sigma_rate"], size = size)
    #sigma = 1 / gamma.rvs(hyperparameters["sigma_shape"], scale = hyperparameters["sigma_scale"], size = size)
    parameters = zeros((3, size))
    parameters[0, :] = nu
    parameters[1, :] = xi
    parameters[2, :] = sigma
    return parameters

## Prior hyperparameters
hyperparameters = { \
        #"nu_shape": 1, \
        #"nu_scale": 3, \
        "nu_rate": 0.2, \
        #"xi_mean": 0, \
        #"xi_sd": 10, \
        "xi_rate": 0.5, \
        #"sigma_shape": 1, \
        #"sigma_scale": 3, \
        "sigma_rate": 0.2}

modeltheta = ParameterModel(name = "Athletics records", dimension = 3)
modeltheta.setHyperparameters(hyperparameters)
modeltheta.setPriorlogdensity(logdprior)
modeltheta.setPriorgenerator(rprior)
modeltheta.setParameterNames(["expression(nu)", "expression(xi)", "expression(sigma)"])
modeltheta.setTransformation(["log", "log", "log"])
modeltheta.plottingInstructions = ["No observations", "No evidence"]
modeltheta.additionalPlots = """
filteringDF <- data.frame(observations)
filteringDF$year <- 1976:2010
filteringDF$xfiltered <- filteredfirststate
filteringDF$CIlow <- filteredCIlow
filteringDF$CIhigh <- filteredCIhigh
names(filteringDF) <- c("y1", "y2", "year", "x", "CIlow", "CIhigh")
g <- ggplot(filteringDF, aes(x = year))
g <- g + geom_line(aes(y = x), size = 1, linetype = 2, colour = "blue")
g <- g + geom_line(aes(y = CIlow), size = 1, colour = "blue")
g <- g + geom_line(aes(y = CIhigh), size = 1, colour = "blue")
g <- g + geom_point(aes(y = y1), size = 4, colour = "green")
g <- g + geom_point(aes(y = y2), size = 4, colour = "red")
g <- g + ylab("Time (seconds)") + ylim(450, 540)
print(g)
smoothingDF <- data.frame(observations)
smoothingDF$year <- 1976:2010
smoothingDF$xsmoothed <- smoothedmeansfirststate35
smoothingDF$CIlow <- smoothedmeansCIlow35
smoothingDF$CIhigh <- smoothedmeansCIhigh35
names(smoothingDF) <- c("y1", "y2", "year", "x", "CIlow", "CIhigh")
g <- ggplot(smoothingDF, aes(x = year))
g <- g + geom_line(aes(y = x), size = 1, linetype = 2, colour = "black")
g <- g + geom_line(aes(y = CIlow), size = 1, colour = "blue")
g <- g + geom_line(aes(y = CIhigh), size = 1, colour = "blue")
g <- g + geom_point(aes(y = y1), size = 4, colour = "green")
g <- g + geom_point(aes(y = y2), size = 4, colour = "red")
g <- g + ylab("Time (seconds)")
g <- g + ylim(450, 540)
print(g)
"""

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
Rpriornu = "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["nu_rate"]
Rpriorxi = "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["xi_rate"]
Rpriorsigma = "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["sigma_rate"]
#Rpriorsigma = InverseGammaTemplate % (hyperparameters["sigma_shape"], hyperparameters["sigma_scale"])
modeltheta.setRprior([Rpriornu, Rpriorxi, Rpriorsigma])

