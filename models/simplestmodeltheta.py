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
        ones, mean, average, prod, log, array
from scipy.stats import norm, truncnorm, gamma
from src.models import ParameterModel


#### See src/models.py for explanations about the model functions.
#### most functions here take transformed parameters, but rprior returns
#### untransformed parameters

def safelogdlogit(x):
    #y = x.copy()
    indTooBig = (x > 10)
    indTooSmall = (x < -10)
    indProblem = indTooBig | indTooSmall
#    y[1 - indProblem] = log(1 + exp(x[1 - indProblem]))
    result = x - 2 * log(1 + exp(x))
    result[indProblem] = -10**10
    return result
def logdprior(parameters, hyperparameters):
    """ Takes transformed parameters.  When the parameter is transformed, 
    a jacobian appears in the formula.
    """
    # the following is the log density of Y = logit(U) when U is Uniform(0,1)
    rho_part = safelogdlogit(array([parameters[0]]))
    #parameters[0] - 2 * log(1 + exp(parameters[0]))
    return rho_part

def rprior(size, hyperparameters):
    """ returns untransformed parameters """
    rho = random.uniform(size = size, low = 0.01, high = 0.99) 
    #sigma = 1 / gamma.rvs(hyperparameters["sigma_shape"], scale = hyperparameters["sigma_scale"], size = size)
    #tau = 1 / gamma.rvs(hyperparameters["tau_shape"], scale = hyperparameters["tau_scale"], size = size)
    parameters = zeros((1, size))
    parameters[0, :] = rho
    return parameters

## Prior hyperparameters
hyperparameters = {}

modeltheta = ParameterModel(name = "Simplest model theta", dimension = 1)
modeltheta.setHyperparameters(hyperparameters)
modeltheta.setPriorlogdensity(logdprior)
modeltheta.setPriorgenerator(rprior)
modeltheta.setParameterNames(["expression(rho)"])
modeltheta.setTransformation(["logit"])
modeltheta.setRtruevalues([0.8])
modeltheta.additionalPlots = """
if ("predictedlowquantile" %in% ls() && "predictedhiquantile" %in% ls()){
    g <- qplot(x = 1:T, y = observations, geom = "line")
    g <- g + geom_line(aes(y = predictedlowquantile), colour = "red")
    g <- g + geom_line(aes(y = predictedhiquantile), colour = "green")
    g <- g + xlab("time") + ylab("observations")
    print(g)
}
if (exists("truestates")){
    if (length(truestates) > 1){
      g <- qplot(x = 1:T, y = truestates, geom = "line")
      g <- g + geom_line(aes(y = filteredhiddenstate), colour = "red")
      g <- g + xlab("time") + ylab("true states")
      print(g)
    }
}
"""
uniformprior = \
"""
priorfunction <- function(x){
    return(1)
}
"""
modeltheta.setRprior([uniformprior])

#InverseGammaTemplate % (hyperparameters["tau_shape"], hyperparameters["tau_scale"])])
#modeltheta.additionalPlots = """
#if ("predictedlowquantile" %in% ls() && "predictedhiquantile" %in% ls()){
#    g <- qplot(x = 1:T, y = observations, geom = "line")
#    g <- g + geom_line(aes(y = predictedlowquantile), colour = "red")
#    g <- g + geom_line(aes(y = predictedhiquantile), colour = "green")
#    g <- g + xlab("time") + ylab("observations")
#    print(g)
#}
#"""
#

