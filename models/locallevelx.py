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
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, newaxis, \
        argsort, var, cumsum, searchsorted
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
from src.models import SSM

################################################################
# Linear Gaussian Model
# X = X + sigma Epsilon
# Y = X + tau Eta
# X_0 ~ N(0, 1)
# parameters[0, :] = sigma ^ 2
# parameters[1, :] = tau ^ 2
################################################################

### these functions take untransformed parameters as arguments

#### See src/models.py for explanations about the model functions.
def firstStateGenerator(parameters, size):
    return random.normal(size = size, loc = 0, scale = 1)[:, newaxis]
def observationGenerator(states, parameters):
    return random.normal(size = states.shape[0], loc = states[:, 0], scale = sqrt(parameters[1]))[:, newaxis]
def transitionAndWeight(states, y, parameters, t):
    code = \
    """
    float tempmeasure1;
    float tempmeasure2;
    float temptransition;
    for (int j = 0; j < Ntheta; j++)
    {
        tempmeasure1 = -0.9189385 - 0.5 * log(parameters(1, j));
        tempmeasure2 = -0.5 / parameters(1, j);
        temptransition = sqrt(parameters(0, j));
        for(int k = 0; k < Nx; k++)
        {
            states(k, 0, j) = states(k, 0, j) + temptransition * noise(k, j);
            weights(k, j) = tempmeasure1 + 
            tempmeasure2 * ((double) y(0) - states(k, 0, j)) * ((double) y(0) - states(k, 0, j));
        }
    }
    """
    y = array([y])
    Nx = states.shape[0]
    Ntheta = states.shape[2]
    weights = zeros((Nx, Ntheta))
    noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1)
    weave.inline(code,['Nx', 'Ntheta', 'states', 'y', 'parameters', 'noise', 'weights'], type_converters=weave.converters.blitz, libraries = ["m"])
    return {"states": states , "weights": weights}

modelx = SSM(name = "Linear Gaussian model x", xdimension = 1, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.parameters = array([0.5, 0.1])
def firststate(xparticles, thetaparticles, t):
    return xparticles[:, 0, :]
modelx.setFiltering({"firststate": firststate})

modelx.setRLinearGaussian(\
"""
dlm <- list("FF" = 1, "GG" = 1, "V" = %.3f, "W" = %.3f,
             "m0" = 0, "C0" = 1)
""" % (modelx.parameters[1], modelx.parameters[0]))
def predictionSquaredObservations(xparticles, thetaweights, thetaparticles, t):
    Nx = xparticles.shape[0]
    Ntheta = xparticles.shape[2]
    result = zeros(3)
    observations = zeros(Nx * Ntheta)
    weightobs = zeros(Nx * Ntheta)
    for j in range(Ntheta):
        observations[(Nx * j):(Nx * (j+1))] = \
                observationGenerator(xparticles[..., j], thetaparticles[:, j]).reshape(Nx)
        weightobs[(Nx * j):(Nx * (j+1))] = repeat(thetaweights[j], repeats = Nx)
    observations = power(observations, 2)
    weightobs = weightobs / sum(weightobs)
    obsmean = average(observations, weights = weightobs)
    ind = argsort(observations)
    observations = observations[ind]
    weightobs = weightobs[ind]
    cumweightobs = cumsum(weightobs)
    quantile5 = observations[searchsorted(cumweightobs, 0.05)]
    quantile95 = observations[searchsorted(cumweightobs, 0.95)]
    result[0] = obsmean
    result[1] = quantile5
    result[2] = quantile95
    return result

def predictionHiddenstate(xparticles, thetaweights, thetaparticles, t):
    Nx = xparticles.shape[0]
    Ntheta = xparticles.shape[2]
    result = zeros(3)
    predictedstate = zeros(Nx * Ntheta)
    weight = zeros(Nx * Ntheta)
    for j in range(Ntheta):
        predictedstate[(Nx * j):(Nx * (j+1))] = xparticles[..., 0, j]
        weight[(Nx * j):(Nx * (j+1))] = repeat(thetaweights[j], repeats = Nx)
    weight = weight / sum(weight)
    xmean = average(predictedstate, weights = weight)
    ind = argsort(predictedstate)
    predictedstate = predictedstate[ind]
    weight = weight[ind]
    cumweight = cumsum(weight)
    quantile5 = predictedstate[searchsorted(cumweight, 0.05)]
    quantile95 = predictedstate[searchsorted(cumweight, 0.95)]
    result[0] = xmean
    result[1] = quantile5
    result[2] = quantile95
    return result
def predictionObservations(xparticles, thetaweights, thetaparticles, t):
    Nx = xparticles.shape[0]
    Ntheta = xparticles.shape[2]
    result = zeros(3)
    observations = zeros(Nx * Ntheta)
    weightobs = zeros(Nx * Ntheta)
    for j in range(Ntheta):
        observations[(Nx * j):(Nx * (j+1))] = \
                observationGenerator(xparticles[..., j], thetaparticles[:, j]).reshape(Nx)
        weightobs[(Nx * j):(Nx * (j+1))] = repeat(thetaweights[j], repeats = Nx)
    weightobs = weightobs / sum(weightobs)
    obsmean = average(observations, weights = weightobs)
    ind = argsort(observations)
    observations = observations[ind]
    weightobs = weightobs[ind]
    cumweightobs = cumsum(weightobs)
    quantile5 = observations[searchsorted(cumweightobs, 0.05)]
    quantile95 = observations[searchsorted(cumweightobs, 0.95)]
    result[0] = obsmean
    result[1] = quantile5
    result[2] = quantile95
    return result

modelx.setPrediction([{"function": predictionSquaredObservations, "dimension": 3, "name": "squaredobs"}, \
        {"function": predictionHiddenstate, "dimension": 3, "name": "hiddenstate"}, \
        {"function": predictionObservations, "dimension": 3, "name": "obs"}])


Rmarginals = \
"""
temptrueloglikelihood <- function(theta){
    somedlm <- dlm
    somedlm["V"] <- theta[2]
    somedlm["W"] <- theta[1]
    return(KFLL(observations, somedlm))
}
trueloglikelihood <- function(theta1, theta2){
  theta <- cbind(theta1, theta2)
  return(apply(X=theta, MARGIN=1, FUN=temptrueloglikelihood))
}
priorfunction <- function(x){
    shape <- 1.00000 
    scale <- 1.00000
    return(scale**shape / gamma(shape) * x**(- shape - 1) * exp(-scale / x))
}
truelogposterior <- function(theta1, theta2){
  prior <- log(priorfunction(theta1)) + log(priorfunction(theta2))
  theta <- cbind(theta1, theta2)
  return(prior + apply(X=theta, MARGIN=1, FUN=temptrueloglikelihood))
}
range1 <- range(thetahistory[1,1,])
range2 <- range(thetahistory[1,2,])
x = seq(from = range1[1], to = range1[2], length.out = 6)
y = seq(from = range2[1], to = range2[2], length.out = 6)
z = outer(X=x,Y=y,FUN=truelogposterior)
maxlogposterior <- max(z)
tmpunnormmarginal1 <- function(theta1){
  temp <- function(theta2){
    exp(truelogposterior(theta1, theta2) - maxlogposterior)
  }
  integrate(temp, lower = range2[1], upper = range2[2])$value
}
unnormmarginal1 <- function(theta1) sapply(X=theta1, FUN=tmpunnormmarginal1)
normconstant1 <- integrate(unnormmarginal1, lower = range1[1], upper = range1[2])$value
normmarginal1 <- function(theta1) unnormmarginal1(theta1) / normconstant1
tmpunnormmarginal2 <- function(theta2){
  temp <- function(theta1){
    exp(truelogposterior(theta1, theta2) - maxlogposterior)
  }
  integrate(temp, lower = range1[1], upper = range1[2])$value
}
unnormmarginal2 <- function(theta2) sapply(X=theta2, FUN=tmpunnormmarginal2)
normconstant2 <- integrate(unnormmarginal2, lower = range2[1], upper = range2[2])$value
normmarginal2 <- function(theta2) unnormmarginal2(theta2) / normconstant2
marginals <- c(normmarginal1, normmarginal2)
"""
modelx.setRmarginals(Rmarginals)

