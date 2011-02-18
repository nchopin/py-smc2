#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from numpy import random, power, sqrt, exp, zeros, zeros_like, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32
from scipy.stats import norm, truncnorm, gamma
from src.models import ParameterModel


#def invgamma_logpdf(logx, shape, scale):
#    return (-shape - 1) * logx - scale / exp(logx)

def logdprior(parameters, hyperparameters):
    """ Takes transformed parameters.  When the parameter is transformed, 
    a jacobian appears in the formula.
    """
    mu_part = - 0.5 / (hyperparameters["mu_sd"]**2) * (parameters[0] - hyperparameters["mu_mean"])**2
    beta_part = - 0.5 / (hyperparameters["beta_sd"]**2) * (parameters[1] - hyperparameters["beta_mean"])**2
    xi_part = parameters[2] - hyperparameters["xi_rate"] * exp(parameters[2])
    omega2_part = parameters[3] - hyperparameters["omega2_rate"] * exp(parameters[3])
    lamb_part = parameters[4] - hyperparameters["lambda_rate"] * exp(parameters[4])
    return mu_part + beta_part + xi_part + omega2_part + lamb_part


def rprior(size, hyperparameters):
    """ returns untransformed parameters """
    mu = norm.rvs(size = size, loc = hyperparameters["mu_mean"], scale = hyperparameters["mu_sd"])
    beta = norm.rvs(size = size, loc = hyperparameters["beta_mean"], scale = hyperparameters["beta_sd"])
    xi = random.exponential(scale = 1 / hyperparameters["xi_rate"], size = size)
    omega2 = random.exponential(scale = 1 / hyperparameters["omega2_rate"], size = size)
    lamb = random.exponential(scale = 1 / hyperparameters["lambda_rate"], size = size)
    parameters = zeros((5, size))
    parameters[0, :] = mu
    parameters[1, :] = beta
    parameters[2, :] = xi
    parameters[3, :] = omega2
    parameters[4, :] = lamb
    return parameters

def rInitDistribution(size):
    """ returns untransformed parameters """
    mu = norm.rvs(size = size, loc = 0, scale = sqrt(2))
    beta = norm.rvs(size = size, loc = 0, scale = sqrt(2))
    xi = random.exponential(scale = 1, size = size)
    omega2 = random.exponential(scale = 1, size = size)
    lamb = random.exponential(scale = 1, size = size)
    parameters = zeros((5, size))
    parameters[0, :] = mu
    parameters[1, :] = beta
    parameters[2, :] = xi
    parameters[3, :] = omega2
    parameters[4, :] = lamb
    return parameters

def dInitDistribution(parameters):
    """ takes transformed parameters """
    mu_part = - 0.5 / (2) * (parameters[0] - 0)**2
    beta_part = - 0.5 / (2) * (parameters[1] - 0)**2
    xi_part = parameters[2] - 1 * exp(parameters[2])
    omega2_part = parameters[3] - 1 * exp(parameters[3])
    lamb_part = parameters[4] - 1 * exp(parameters[4])
    return mu_part + beta_part + xi_part + omega2_part + lamb_part


hyperparameters = { \
        "mu_mean": 0, "mu_sd": sqrt(2), \
        "beta_mean": 0, "beta_sd": sqrt(2), \
        "xi_rate": 0.2, "omega2_rate": 0.2, \
        "lambda_rate": 1}


modeltheta = ParameterModel(name = "SV one-factor", dimension = 5)
modeltheta.setHyperparameters(hyperparameters)
modeltheta.setPriorlogdensity(logdprior)
modeltheta.setPriorgenerator(rprior)
modeltheta.setInitDistribution(rInitDistribution, dInitDistribution)
modeltheta.setParameterNames(["expression(mu)", "expression(beta)", \
        "expression(xi)", "expression(omega^2)", "expression(lambda)"])
modeltheta.setTransformation(["none", "none", "log", "log", "log"])
modeltheta.priorandtruevaluesspecified = False
modeltheta.Rparameterstruevalues = "trueparameters <- c(0, 0, 0.5, 0.0625, 0.01)"
modeltheta.setRprior(["priorfunction <- function(x) dnorm(x, sd = %.5f)" % hyperparameters["mu_sd"], \
                            "priorfunction <- function(x) dnorm(x, sd = %.5f)" % hyperparameters["beta_sd"], \
                            "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["xi_rate"], \
                            "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["omega2_rate"], \
                            "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["lambda_rate"]])








