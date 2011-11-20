#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from numpy import random, power, sqrt, exp, zeros, zeros_like, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32
from scipy.stats import norm, truncnorm, gamma
from src.parametermodel import ParameterModel


#def invgamma_logpdf(logx, shape, scale):
#    return (-shape - 1) * logx - scale / exp(logx)

def logdprior(parameters, hyperparameters):
    """ takes transformed parameters """
    xi_part = parameters[0] - hyperparameters["xi_rate"] * exp(parameters[0])
    omega2_part = parameters[1] - hyperparameters["omega2_rate"] * exp(parameters[1])
    lamb_part = parameters[2] - hyperparameters["lambda_rate"] * exp(parameters[2])
    return xi_part + omega2_part + lamb_part


def rprior(size, hyperparameters):
    """ returns untransformed parameters """
    xi = random.exponential(scale = 1 / hyperparameters["xi_rate"], size = size)
    omega2 = random.exponential(scale = 1 / hyperparameters["omega2_rate"], size = size)
    lamb = random.exponential(scale = 1 / hyperparameters["lambda_rate"], size = size)
    parameters = zeros((3, size))
    parameters[0, :] = xi
    parameters[1, :] = omega2
    parameters[2, :] = lamb
    return parameters

hyperparameters = { \
        "xi_rate": 0.2, "omega2_rate": 0.2, \
        "lambda_rate": 1}

def rInitDistribution(size):
    """ returns untransformed parameters """
    xi = random.exponential(scale = 1, size = size)
    omega2 = random.exponential(scale = 1, size = size)
    lamb = random.exponential(scale = 1, size = size)
    parameters = zeros((3, size))
    parameters[0, :] = xi
    parameters[1, :] = omega2
    parameters[2, :] = lamb
    return parameters

def dInitDistribution(parameters):
    """ takes transformed parameters """
    xi_part = parameters[0] - 1 * exp(parameters[0])
    omega2_part = parameters[1] - 1 * exp(parameters[1])
    lamb_part = parameters[2] - 1 * exp(parameters[2])
    return xi_part + omega2_part + lamb_part




modeltheta = ParameterModel(name = "SV one-factor", dimension = 3)
modeltheta.setHyperparameters(hyperparameters)
modeltheta.setPriorlogdensity(logdprior)
modeltheta.setPriorgenerator(rprior)
modeltheta.setInitDistribution(rInitDistribution, dInitDistribution)
modeltheta.setParameterNames(["expression(xi)", "expression(omega^2)", "expression(lambda)"])
modeltheta.setTransformation(["log", "log", "log"])
modeltheta.priorandtruevaluesspecified = True
modeltheta.Rparameterstruevalues = "trueparameters <- c(0.5, 0.0625, 0.01)"
modeltheta.Rpriorfunctions = ["priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["xi_rate"], \
                            "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["omega2_rate"], \
                            "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["lambda_rate"]]





