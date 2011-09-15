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
# Periodic Gaussian Model
# X = 0.5 X + 25 X / (1 + X^2) + 8 cos(1.2(t - 1)) + sigma_w W_t
# Y = X^2 / 20 + sigma_v V_t
# X_0 ~ N(0,2)
# parameters[0, :] = sigma_w ^ 2
# parameters[1, :] = sigma_v ^ 2
################################################################

### these functions take untransformed parameters as arguments

#### See src/models.py for explanations about the model functions.
def firstStateGenerator(parameters, size):
    #return random.normal(size = size, loc = 0, scale = sqrt(2))[:, newaxis]
    return zeros((size, 1)) + 0.1
def observationGenerator(states, parameters):
    return random.normal(size = states.shape[0], loc = power(states[:, 0], 2) / 20, \
                    scale = sqrt(parameters[1]))[:, newaxis]
def transitionAndWeight(states, y, parameters, t):
    code = \
    """
    float tempmeasure1;
    float tempmeasure2;
    float temptransition;
    float term;
    for (int j = 0; j < Ntheta; j++)
    {
        tempmeasure1 = -0.9189385 - 0.5 * log(parameters(1, j));
        tempmeasure2 = -0.5 / parameters(1, j);
        temptransition = sqrt(parameters(0, j));
        for(int k = 0; k < Nx; k++)
        {
            states(k, 0, j) = 0.5 * states(k, 0, j) + 25 * states(k, 0, j) /
            (1 + states(k, 0, j) * states(k, 0, j)) + 8 * cos(1.2 * (t(0) - 1))
            + temptransition * noise(k, j);
            term = (double) y(0) - (states(k, 0, j) * states(k, 0, j) / 20);
            weights(k, j) = tempmeasure1 + tempmeasure2 * term * term;
        }
    }
    """
    y = array([y])
    t = array([t])
    Nx = states.shape[0]
    Ntheta = states.shape[2]
    weights = zeros((Nx, Ntheta))
    noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1)
    weave.inline(code,['Nx', 'Ntheta', 'states', 'y', 't', 'parameters', 'noise', 'weights'], type_converters=weave.converters.blitz, libraries = ["m"])
    return {"states": states , "weights": weights}

modelx = SSM(name = "Periodic Gaussian model x", xdimension = 1, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.parameters = array([10, 1])




