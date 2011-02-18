#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, zeros_like, newaxis
from numpy import sum as numpysum
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
from src.models import SSM
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
# Here mu and beta are constants = 0.
# parameters[0, :] = xi
# parameters[1, :] = omega^2
# parameters[2, :] = lambda
############################################

### these functions take untransformed parameters as arguments

def firstStateGenerator(parameters, size):
    first_state = zeros((size, 2))
    first_state[:, 1] = random.gamma(size = size, shape = (parameters[0]**2) / parameters[1], scale = (parameters[1] / parameters[0]))
    return first_state
def observationGenerator(states, parameters):
    return random.normal(size = states.shape[0], loc = 0, scale = sqrt(abs(states[:, 0])))[:, newaxis]

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
    exponentialscale = parameters(1) / parameters(0);
    for(int indexk = 0; indexk < Nx; indexk ++){
        sum1 = 0.;
        sum2 = 0.;
        k = allK(indexk);
        for (int indexl = 0; indexl < k; indexl ++){
            auxiliaryE = - log(alluniforms(currentk)) * exponentialscale;
            auxiliaryC = alluniforms(currentk + 1);
            currentk += 2;
            sum1 = sum1 + auxiliaryE;
            sum2 = sum2 + auxiliaryE * exp(-parameters(2) * auxiliaryC);
        }
        zs = exp(-parameters(2)) * states(indexk, 1) + sum2;
        vs = (1 / parameters(2)) * (states(indexk, 1) - zs + sum1);
        states(indexk, 0) = vs;
        states(indexk, 1) = zs;
        tempmeasure1 = y(0);
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
    # --------------
    poissonparameters = parameters[2, :] * (parameters[0, :]**2) / parameters[1, :]
    poissonparameters = repeat(poissonparameters[:,newaxis], Nx, axis = 1)
    for indextheta in range(Ntheta):
        allK = random.poisson(lam = poissonparameters[indextheta,:])
        allK = array(allK).reshape(Nx)
        sumK = sum(allK)
        alluniforms = random.uniform(size = 2 * sumK)
        subresults = subtransitionAndWeight(states[..., indextheta], y, parameters[:, indextheta], \
                         alluniforms, allK)
        newstates[..., indextheta] = subresults["states"]
        weights[..., indextheta] = subresults["weights"]
    # --------------
    return {"states": newstates , "weights": weights}



modelx = SSM("SV one-factor", xdimension = 2, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.parameters = array([0.5, 0.0625, 0.01])





