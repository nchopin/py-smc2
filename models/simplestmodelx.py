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
        array, float32, int32, zeros_like, newaxis
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
from src.models import SSM

################################################################
# Simplest Linear Gaussian Model
# X = rho * X + sigma Epsilon
# Y = X + tau Eta
# X_0 ~ N(0, 1)
# parameters[0, :] = rho
# sigma = 0.5
# tau = 0.5
################################################################

SIGMA = 0.5
TAU = 0.5
SIGMA2 = SIGMA * SIGMA
TAU2 = TAU * TAU
rho = 0.8

### these functions take untransformed parameters as arguments
#### See src/models.py for explanations about the model functions.
def firstStateGenerator(parameters, size):
    return random.normal(size = size, loc = 0, scale = 1)[:, newaxis]
def observationGenerator(states, parameters):
    return random.normal(size = states.shape[0], loc = states[:, 0], scale = TAU)[:, newaxis]
def transitionAndWeight(states, y, parameters, t):
    code = \
    """
    float tempmeasure1;
    float tempmeasure2;
    float temptransition;
    for (int j = 0; j < Ntheta; j++)
    {
        tempmeasure1 = -0.9189385 - 0.5 * log(%(TAU2)s);
        tempmeasure2 = -0.5 / (%(TAU2)s);
        temptransition = %(SIGMA)s;
        for(int k = 0; k < Nx; k++)
        {
            states(k, 0, j) = parameters(0, j) * states(k, 0, j) + temptransition * noise(k, j);
            weights(k, j) = tempmeasure1 + 
            tempmeasure2 * ((double) y(0) - states(k, 0, j)) * ((double) y(0) - states(k, 0, j));
        }
    }
    """ % {"TAU2": TAU2, "SIGMA": SIGMA}
    y = array([y])
    Nx = states.shape[0]
    Ntheta = states.shape[2]
    weights = zeros((Nx, Ntheta))
    noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1)
    weave.inline(code,['Nx', 'Ntheta', 'states', 'y', 'parameters', 'noise', 'weights'], \
            type_converters=weave.converters.blitz, libraries = ["m"])
    return {"states": states , "weights": weights}

modelx = SSM(name = "Simplest model x", xdimension = 1, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.parameters = array([rho])

def firststate(xparticles, thetaparticles, t):
    return xparticles[:, 0, :]
modelx.setFiltering({"firststate": firststate})
# this model is a linear gaussian model, we can specify it
modelx.setRLinearGaussian(\
"""
dlm <- list("FF" = 1, "GG" = %.3f, "V" = %.3f, "W" = %.3f,
             "m0" = 0, "C0" = 1)
""" % (modelx.parameters[0], SIGMA2, TAU2))

#Rtruelikelihood = \
#"""
#temptrueloglikelihood <- function(theta){
#    somedlm <- dlm
#    somedlm["GG"] <- theta
#    return(KFLL(observations, somedlm))
#}
#trueloglikelihood <- function(theta){
#    return(sapply(X= theta, FUN= temptrueloglikelihood))
#}
#trueunnormlikelihood <- function(theta) exp(trueloglikelihood(theta))
#normlikelihood <- integrate(f = trueunnormlikelihood, lower = 0, upper = 1)$value
#truelikelihood <- function(theta) trueunnormlikelihood(theta) / normlikelihood
#"""
#modelx.setRmarginals([Rtruelikelihood])




