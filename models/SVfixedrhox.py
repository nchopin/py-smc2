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
        array, zeros_like, newaxis
from numpy import sum as numpysum
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
from src.models import SSM
from snippets.localfolder import get_path

############################################
## Stochastic volatility : multi-factor model with fixed rho
# 
# Y_t = mu + beta * v_t + v_t^(0.5) * epsilon_t + ...
# X_t = (v_t, v(1)_t, z(1)_t, v(2)_t, z(2)_t)
# v_t = v(1)_t + v(2)_t
# Each v(i)_t, z(i)_t:
## v_t+1 = lambda^(-1) ( z_t - z_{t+1} + sum_j=1^k e_j )
## z_t+1 = e^(-lambda) * z_t + sum_j=1^k(1 - exp(-lambda(t + 1 - c_j))) e_j
## k ~ Poisson(lambda * xi^2 / omega^2)
## c_{1:k} ~ Uniform(t, t + 1)
## e_{1:k} ~ Exp(xi / omega^2) (rate parameter)
#
# v_0 does not matter
# z(i)_0 ~ Gamma(w(i) * xi^2 / omega^2, xi / omega^2)
#
# parameters[0, :] = mu
# parameters[1, :] = beta
# parameters[2, :] = xi
# parameters[3, :] = omega^2
# parameters[4, :] = lambda(1)
# parameters[5, :] = lambda(2) - lambda(1)
# parameters[6, :] = w(1)
#
# states[:, 0] = v
# states[:, 1] = v(1)
# states[:, 2] = z(1)
# states[:, 3] = v(2)
# states[:, 4] = z(2)
# states[:, 5] = sum of e_{1, j}
# states[:, 6] = sum of e_{2, j}
############################################

### these functions take untransformed parameters as arguments

def firstStateGenerator(parameters, size):
    first_state = zeros((size, 5))
    first_state[:, 2] = random.gamma(size = size, shape = parameters[6] * (parameters[2]**2) / parameters[3], scale = (parameters[3] / parameters[2]))
    try:
        first_state[:, 4] = random.gamma(size = size, shape = (1 - parameters[6]) * (parameters[2]**2) / parameters[3], scale = (parameters[3] / parameters[2]))
    except:
        print "pb generating gamma"
        print "parameters"
        print parameters
        #raw_input("go on?")
        raise ValueError("ERROR: generating the initial gamma variables") 
    return first_state
def observationGenerator(states, parameters):
## not implemented
    return random.normal(size = states.shape[0], loc = 0, scale = 1)

def subtransitionAndWeight(states, y, parameters, alluniforms1, allK1, alluniforms2, allK2):
    code = \
    """
    float tempmeasure1 = 0.;
    float zs1 = 0.;
    float vs1 = 0.;
    float zs2 = 0.;
    float vs2 = 0.;
    int currentk1 = 0;
    int currentk2 = 0;
    float exponentialscale = 0.;
    float sum11 = 0.;
    float sum21 = 0.;
    float sum12 = 0.;
    float sum22 = 0.;
    int k1 = 0;
    int k2 = 0;
    float auxiliaryE1 = 0.;
    float auxiliaryC1 = 0.;
    float auxiliaryE2 = 0.;
    float auxiliaryC2 = 0.;
    exponentialscale = parameters(3) / parameters(2);
    for(int indexk = 0; indexk < Nx; indexk ++){
        sum11 = 0.;
        sum21 = 0.;
        k1 = allK1(indexk);
        for (int indexl = 0; indexl < k1; indexl ++){
            auxiliaryE1 = - log(alluniforms1(currentk1)) * exponentialscale;
            auxiliaryC1 = alluniforms1(currentk1 + 1);
            currentk1 += 2;
            sum11 = sum11 + auxiliaryE1;
            sum21 = sum21 + auxiliaryE1 * exp(-parameters(4) * auxiliaryC1);
        }
        sum12 = 0.;
        sum22 = 0.;
        k2 = allK2(indexk);
        for (int indexl = 0; indexl < k2; indexl ++){
            auxiliaryE2 = - log(alluniforms2(currentk2)) * exponentialscale;
            auxiliaryC2 = alluniforms2(currentk2 + 1);
            currentk2 += 2;
            sum12 = sum12 + auxiliaryE2;
            sum22 = sum22 + auxiliaryE2 * exp(-(parameters(4) + parameters(5)) * auxiliaryC2);
        }
        zs1 = exp(-parameters(4)) * states(indexk, 2) + sum21;
        vs1 = (1 / parameters(4)) * (states(indexk, 2) - zs1 + sum11);
        states(indexk, 1) = vs1;
        states(indexk, 2) = zs1;
        zs2 = exp(- (parameters(4) + parameters(5))) * states(indexk, 4) + sum22;
        vs2 = (1 / (parameters(4) + parameters(5))) * (states(indexk, 4) - zs2 + sum12);
        states(indexk, 3) = vs2;
        states(indexk, 4) = zs2;
        states(indexk, 0) = vs1 + vs2;
        tempmeasure1 = y(0) - (parameters(0) + parameters(1) * states(indexk, 0));
        weights(indexk) = -0.9189385 - 0.5 * log(states(indexk, 0)) - 0.5 / states(indexk, 0) *
                (tempmeasure1 * tempmeasure1);
    }
    """
    y = array([y])
    Nx = states.shape[0]
    weights = zeros(Nx)
    weave.inline(code,['Nx', 'states', 'y', 'parameters', 'weights', 'alluniforms1', 'allK1', 'alluniforms2', 'allK2'], \
            type_converters=weave.converters.blitz, libraries = ["m"])
    return {"states": states , "weights": weights}


def transitionAndWeight(states, y, parameters, t):
    Nx = states.shape[0]
    Ntheta = states.shape[2]
    weights = zeros((Nx, Ntheta))
    newstates = zeros_like(states)
    # --------------
    poissonparameters1 = parameters[6, :] * parameters[4, :] * (parameters[2, :]**2) / parameters[3, :]
    poissonparameters2 = (1 - parameters[6, :]) * (parameters[4, :] + parameters[5, :]) * (parameters[2, :]**2) / parameters[3, :]
    poissonparameters1 = repeat(poissonparameters1[:,newaxis], Nx, axis = 1)
    poissonparameters2 = repeat(poissonparameters2[:,newaxis], Nx, axis = 1)
    for indextheta in range(Ntheta):
        allK1 = random.poisson(lam = poissonparameters1[indextheta,:])
        allK1 = array(allK1).reshape(Nx)
        sumK1 = sum(allK1)
        allK2 = random.poisson(lam = poissonparameters2[indextheta,:])
        allK2 = array(allK2).reshape(Nx)
        sumK2 = sum(allK2)
        alluniforms1 = random.uniform(size = 2 * sumK1)
        alluniforms2 = random.uniform(size = 2 * sumK2)
        subresults = subtransitionAndWeight(states[..., indextheta], y, parameters[:, indextheta], \
                         alluniforms1, allK1, alluniforms2, allK2)
        newstates[..., indextheta] = subresults["states"]
        weights[..., indextheta] = subresults["weights"]
    # --------------
    return {"states": newstates , "weights": weights}

modelx = SSM("SV multi-factor", xdimension = 5, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
# stupid values in this case, we don't simulate datasets from this model anyway
modelx.parameters = array([0, 0, 0., 0., 0., 0, 0.])





