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
        array, zeros_like, newaxis, float32, int32 
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
from src.models import SSM

############################################
## Population Dynamic model, parameterization of Polansky et al. (2009)
# n_t+1 = n_t + r(1 - (exp(n_t) / K)^theta) + sigma_epsilon Epsilon
# Y_t ~ normal(exp(n_t), sigma_w^2)
# parameters[0, :] = sigma_epsilon ^ 2
# parameters[1, :] = sigma_w ^ 2
# parameters[2, :] = n0
# parameters[3, :] = r
# parameters[4, :] = K
# parameters[5, :] = theta
############################################

### these functions take untransformed parameters as arguments

def firstStateGenerator(parameters, size):
    return repeat(parameters[2], size)[:, newaxis]
def observationGenerator(states, parameters):
    return random.normal(size = states.shape[0], loc = exp(states[:,0]), scale = sqrt(parameters[1]))[:, newaxis]
def transitionAndWeight(states, y, parameters, t):
    code = \
    """
    float tempmeasure1;
    float temptrans;
    for (int j = 0; j < Ntheta; j++){
        tempmeasure1 = -0.9189385 - 0.5 * log(parameters(1, j));
        temptrans = exp(-parameters(5,j) * log(parameters(4, j)));
        for(int k = 0; k < Nx; k++){
                states(k, 0, j) = states(k, 0, j) + 
                parameters(3, j) * (1 - temptrans * exp(parameters(5, j) * states(k, 0, j)))
                + sqrt(parameters(0, j)) * noise(k, j);
                weights(k, j) =  tempmeasure1 - 0.5 / parameters(1, j)
                * ((double) y(0) - exp(states(k, 0, j))) * ((double) y(0) - exp(states(k, 0, j)));
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

modelx = SSM("Population Dynamic model x (M2), reparameterized", xdimension = 1, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.parameters = array([0.05**2, 0.05**2, log(1), 0.18, 1, 2])



