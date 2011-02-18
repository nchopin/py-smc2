#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from numpy import random, power, sqrt, exp, zeros, zeros_like, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32
from scipy.stats import norm, truncnorm, gamma
from src.models import ParameterModel


def logdprior(parameters, hyperparameters):
    """ Takes transformed parameters.  When the parameter is transformed, 
    a jacobian appears in the formula.
    """
    mu_part = - 0.5 / (hyperparameters["mu_sd"]**2) * (parameters[0] - hyperparameters["mu_mean"])**2
    beta_part = - 0.5 / (hyperparameters["beta_sd"]**2) * (parameters[1] - hyperparameters["beta_mean"])**2
    xi_part = parameters[2] - hyperparameters["xi_rate"] * exp(parameters[2])
    omega2_part = parameters[3] - hyperparameters["omega2_rate"] * exp(parameters[3])
    lamb1_part = parameters[4] - hyperparameters["lambda1_rate"] * exp(parameters[4])
    lambadd_part = parameters[5] - hyperparameters["lambdaadd_rate"] * exp(parameters[5])
    w_part = parameters[6] - 2 * log(1 + exp(parameters[6]))
    return mu_part + beta_part + xi_part + omega2_part + lamb1_part + lambadd_part + w_part

def rprior(size, hyperparameters):
    """ returns untransformed parameters """
    mu = norm.rvs(size = size, loc = hyperparameters["mu_mean"], scale = hyperparameters["mu_sd"])
    beta = norm.rvs(size = size, loc = hyperparameters["beta_mean"], scale = hyperparameters["beta_sd"])
    xi = random.exponential(scale = 1 / hyperparameters["xi_rate"], size = size)
    omega2 = random.exponential(scale = 1 / hyperparameters["omega2_rate"], size = size)
    lamb = random.exponential(scale = 1 / hyperparameters["lambda1_rate"], size = size)
    lambadd = random.exponential(scale = 1 / hyperparameters["lambdaadd_rate"], size = size)
    w = random.uniform(size = size)
    parameters = zeros((7, size))
    parameters[0, :] = mu
    parameters[1, :] = beta
    parameters[2, :] = xi
    parameters[3, :] = omega2
    parameters[4, :] = lamb
    parameters[5, :] = lambadd
    parameters[6, :] = w
    return parameters

def rInitDistribution(size):
    """ returns untransformed parameters """
    mu = norm.rvs(size = size, loc = 0, scale = sqrt(2))
    beta = norm.rvs(size = size, loc = 0, scale = sqrt(2))
    xi = random.exponential(scale = 1, size = size)
    omega2 = random.exponential(scale = 1, size = size)
    lamb = random.exponential(scale = 1, size = size)
    lambadd = random.exponential(scale = 2, size = size)
    w = random.uniform(size = size)
    parameters = zeros((7, size))
    parameters[0, :] = mu
    parameters[1, :] = beta
    parameters[2, :] = xi
    parameters[3, :] = omega2
    parameters[4, :] = lamb
    parameters[5, :] = lambadd
    parameters[6, :] = w
    return parameters

def dInitDistribution(parameters):
    """ takes transformed parameters """
    mu_part = - 0.5 / (2) * (parameters[0] - 0)**2
    beta_part = - 0.5 / (2) * (parameters[1] - 0)**2
    xi_part = parameters[2] - 1 * exp(parameters[2])
    omega2_part = parameters[3] - 1 * exp(parameters[3])
    lamb1_part = parameters[4] - 1 * exp(parameters[4])
    lambadd_part = parameters[5] - 0.5 * exp(parameters[5])
    w_part = parameters[6] - 2 * log(1 + exp(parameters[6]))
    return mu_part + beta_part + xi_part + omega2_part + lamb1_part + lambadd_part + w_part


hyperparameters = { \
        "mu_mean": 0, "mu_sd": sqrt(2), \
        "beta_mean": 0, "beta_sd": sqrt(2), \
        "xi_rate": 0.2, "omega2_rate": 0.2, \
        "lambda1_rate": 1, "lambdaadd_rate": 0.5}


modeltheta = ParameterModel(name = "SV multi-factor", dimension = 7)
modeltheta.setHyperparameters(hyperparameters)
modeltheta.setPriorlogdensity(logdprior)
modeltheta.setPriorgenerator(rprior)
modeltheta.setInitDistribution(rInitDistribution, dInitDistribution)
modeltheta.setParameterNames(["expression(mu)", "expression(beta)", \
        "expression(xi)", "expression(omega^2)", "expression(lambda[1])", "expression(lambda[2] - lambda[1])", \
        "expression(w[1])"])
#modeltheta.setTransformation([False, False, True, True, True, True, True, False, False])
modeltheta.setTransformation(["none", "none", "log", "log", "log", "log", "logit"])
modeltheta.additionalPlots = """
observationsDF <- cbind(data.frame(observations ** 2), 1:length(observations))
names(observationsDF) <- c("y", "index")
g <- ggplot(data = observationsDF, aes(x = index, y = y)) 
g <- g + geom_line() + ylab("squared observations")
print(g)
"""
modeltheta.setRprior(["priorfunction <- function(x) dnorm(x, sd = %.5f)" % hyperparameters["mu_sd"], \
                      "priorfunction <- function(x) dnorm(x, sd = %.5f)" % hyperparameters["beta_sd"], \
                      "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["xi_rate"], \
                      "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["omega2_rate"], \
                      "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["lambda1_rate"], \
                      "priorfunction <- function(x) dexp(x, rate = %.5f)" % hyperparameters["lambdaadd_rate"], \
                      "priorfunction <- function(x) 1"])





