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
from numpy import random, power, sqrt, exp, zeros, zeros_like, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, var
from scipy.stats import norm, truncnorm, gamma
from src.models import ParameterModel


def invgamma_logpdf(logx, shape, scale):
    return (-shape - 1) * logx - scale / exp(logx)

def logdprior(parameters, hyperparameters):
    """ Takes transformed parameters.  When the parameter is transformed, 
    a jacobian appears in the formula.
    """
    eps_part = parameters[0] + invgamma_logpdf(parameters[0], hyperparameters["eps_shape"], hyperparameters["eps_scale"])
    w_part = parameters[1] + invgamma_logpdf(parameters[1], hyperparameters["w_shape"], hyperparameters["w_scale"])
    n0_part = - 0.5 / (hyperparameters["n0_sd"]**2) * (exp(parameters[2]) - hyperparameters["n0_mean"])**2
    r_part = parameters[3] - exp(parameters[3]) / hyperparameters["r_scale"]
    K_part = parameters[4] - 0.5 / (hyperparameters["K_sd"]**2) * (exp(parameters[4]) - hyperparameters["K_mean"])**2
    theta_part = parameters[5] - exp(parameters[5]) / hyperparameters["theta_scale"]
    result = eps_part + w_part + n0_part + r_part + K_part + theta_part
    return result


def rprior(size, hyperparameters):
    """ returns untransformed parameters """
    eps = 1 / gamma.rvs(hyperparameters["eps_shape"], scale = hyperparameters["eps_scale"], size = size)
    w = 1 / gamma.rvs(hyperparameters["w_shape"], scale = hyperparameters["w_scale"], size = size)
    n0 = truncnorm.rvs((hyperparameters["n0_a"] - hyperparameters["n0_mean"]) / hyperparameters["n0_sd"], \
            (hyperparameters["n0_b"] - hyperparameters["n0_mean"]) / hyperparameters["n0_sd"], size = size, \
            loc = hyperparameters["n0_mean"], scale = hyperparameters["n0_sd"])
    r = random.exponential(scale = hyperparameters["r_scale"], size = size)
    K = truncnorm.rvs((hyperparameters["K_a"] - hyperparameters["K_mean"]) / hyperparameters["K_sd"], \
            (hyperparameters["K_b"] - hyperparameters["K_mean"]) / hyperparameters["K_sd"], size = size, \
            loc = hyperparameters["K_mean"], scale = hyperparameters["K_sd"])
    theta = random.exponential(scale = hyperparameters["theta_scale"], size = size)
    parameters = zeros((6, size))
    parameters[0, :] = eps
    parameters[1, :] = w
    parameters[2, :] = log(n0)
    parameters[3, :] = r
    parameters[4, :] = K
    parameters[5, :] = theta
    return parameters

## Prior hyperparameters
hyperparameters = { \
        "eps_shape": 2, "eps_scale": 1, \
        "w_shape": 2, "w_scale": 1, \
        "r_scale": 1, \
        "K_mean": 1, "K_sd": 1, "K_a": 0, "K_b": 10**10, \
        "theta_scale": 1, \
        "n0_mean": 1, "n0_sd": 1, "n0_a": 0, "n0_b": 10**10}

def updateHyperparametersWithData(hyperparameters, observations):
    hyperparameters["K_mean"] = mean(observations)
    hyperparameters["K_sd"] = 2 * sqrt(var(observations))
    hyperparameters["n0_mean"] = mean(observations)
    hyperparameters["n0_sd"] = 2 * sqrt(var(observations))
    return hyperparameters


modeltheta = ParameterModel(name = "Population Dynamic model theta (M2), reparameterized", dimension = 6)
modeltheta.setHyperparameters(hyperparameters)
modeltheta.setUpdateHyperparametersWithData(updateHyperparametersWithData)
modeltheta.setPriorlogdensity(logdprior)
modeltheta.setPriorgenerator(rprior)
modeltheta.setParameterNames(["expression(sigma[epsilon]^2)", "expression(sigma[W]^2)", \
        "expression(log(N[0]))", "expression(r)", "expression(K)", "expression(theta)"])
modeltheta.setTransformation(["log", "log", "none", "log", "log", "log"])
modeltheta.setRtruevalues([0.05**2, 0.05**2, log(1), 0.18, 1, 0.2])



