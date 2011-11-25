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
        array, zeros_like, newaxis, \
        argsort, cumsum, searchsorted
from numpy import sum as numpysum
from numpy import max as numpymax
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
from src.statespacemodel import SSM
from snippets.localfolder import get_path

############################################
## Stochastic volatility : one-factor model
# 
# Y_t = mu + beta * v_t + v_t^(0.5) * epsilon_t
# X_t = (v_t, z_t)
# v_t+1 = lambda^(-1) ( z_t - z_{t+1} + sum_j=1^k e_j )
# z_t+1 = e^(-lambda) * z_t + sum_j=1^k(1 - exp(-lambda(t + 1 - c_j))) e_j
# k ~ Poisson(lambda * xi^2 / omega^2)
# c_{1:k} ~ Uniform(t, t + 1)
# e_{1:k} ~ Exp(xi / omega^2) (rate parameter)
#
# v_0 does not matter
# z_0 ~ Gamma(xi^2 / omega^2, xi / omega^2)
#
# parameters[0, :] = mu
# parameters[1, :] = beta
# parameters[2, :] = xi
# parameters[3, :] = omega^2
# parameters[4, :] = lambda
# 0 0 0.5 0.0625 0.01
# lambda * xi^2 / omega^2 = 0.01 * 0.25 / (0.0625) = 0.04
############################################

### these functions take untransformed parameters as arguments

def firstStateGenerator(parameters, size):
    first_state = zeros((size, 2))
    first_state[:, 1] = random.gamma(size = size, shape = (parameters[2]**2) / parameters[3], scale = (parameters[3] / parameters[2]))
    return first_state
def observationGenerator(states, parameters):
    correctedscale = sqrt(abs(states[:, 0])) + 10**(-5) 
    return random.normal(size = states.shape[0], loc = parameters[0] + parameters[1] * states[:, 0], \
            scale = correctedscale)[:, newaxis]
def subtransitionAndWeight(states, y, parameters, alluniforms, allK):
    code = \
    """
    float tempmeasure1 = 0.;
    float zs = 0.;
    float vs = 0.;
    int currentk = 0;
    float exponentialscale = 0.;
    float sum1 = 0.;
    float sum2 = 0.;
    int k = 0;
    float auxiliaryE = 0.;
    float auxiliaryC = 0.;
    exponentialscale = parameters(3) / parameters(2);
    for(int indexk = 0; indexk < Nx; indexk ++){
        sum1 = 0.;
        sum2 = 0.;
        k = allK(indexk);
        for (int indexl = 0; indexl < k; indexl ++){
            auxiliaryE = - log(alluniforms(currentk)) * exponentialscale;
            auxiliaryC = alluniforms(currentk + 1);
            currentk += 2;
            sum1 = sum1 + auxiliaryE;
            sum2 = sum2 + auxiliaryE * exp(-parameters(4) * auxiliaryC);
        }
        zs = exp(-parameters(4)) * states(indexk, 1) + sum2;
        vs = (1 / parameters(4)) * (states(indexk, 1) - zs + sum1);
        states(indexk, 0) = vs;
        states(indexk, 1) = zs;
        tempmeasure1 = y(0) - parameters(0) - parameters(1) * states(indexk, 0);
        weights(indexk) = -0.9189385 - 0.5 * log(states(indexk, 0)) - 0.5 / states(indexk, 0) *
                (tempmeasure1 * tempmeasure1);
    }
    """
    y = array([y])
    Nx = states.shape[0]
    weights = zeros(Nx)
    weave.inline(code,['Nx', 'states', 'y', 'parameters', 'weights', 'alluniforms', 'allK'], \
            type_converters=weave.converters.blitz, libraries = ["m"])
    return {"states": states , "weights": weights}

def transitionAndWeight(states, y, parameters, t):
    Nx = states.shape[0]
    Ntheta = states.shape[2]
    weights = zeros((Nx, Ntheta))
    newstates = zeros_like(states)
    poissonparameters = parameters[4, :] * (parameters[2, :]**2) / parameters[3, :]
    for indextheta in range(Ntheta):
        #print poissonparameters[indextheta,:]
        #allK = random.poisson(lam = poissonparameters[indextheta,:], size = Nx)
        allK = random.poisson(lam = poissonparameters[indextheta], size = Nx)
        #print allK
        #raw_input("")
        #if (numpymax(allK) > 10**4):
        #    print "\n ! number of variables to generate for some theta-particle:", t, indextheta, numpymax(allK)
        # the following prevents computational issues while not changing the results,
        # since parameters such that allK > 10**4 are very very likely to have small weights
        # alternatively this fix can be seen as a modification of the model where
        # the poisson law is replaced by a truncated poisson law
        allK[allK > 10**4] = 10**4
        allK = array(allK).reshape(Nx)
        sumK = sum(allK)
        alluniforms = random.uniform(size = 2 * sumK)
        subresults = subtransitionAndWeight(states[..., indextheta], y, parameters[:, indextheta], \
                         alluniforms, allK)
        newstates[..., indextheta] = subresults["states"]
        weights[..., indextheta] = subresults["weights"]
    return {"states": newstates , "weights": weights}

modelx = SSM("SV one-factor", xdimension = 2, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.setParameters(array([0, 0, 0.5, 0.0625, 0.01]))
modelx.addStateFiltering()
modelx.addStatePrediction()
modelx.addObsPrediction()
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
modelx.addPredictionList([{"function": predictionSquaredObservations, "dimension": 3, "name": "squaredobs"}])





