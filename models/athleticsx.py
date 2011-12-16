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
        array, zeros_like, newaxis, maximum, isnan
from numpy import sum as numpysum
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
from src.statespacemodel import SSM
from snippets.localfolder import get_path

############################################
## Athletics records model
# 
# Y(1:r)_t = ordered_GEV(location = mu_t, scale = sigma, shape = xi)
# X_t = (mu_t, mu'_t)
# mu_t+1 = mu_t + mu'_t + epsilon(1)_t
# mu'_t+1 = mu'_t + epsilon(2)_t
# epsilon = (epsilon(1), epsilon(2)) ~ Normal(0, Q)
# Q = nu^2 * ( 1/3  1/2
#              1/2   1  )
#
# mu_0 ~ Normal(520, 10^2)
#
# parameters[0, :] = nu (and not nu^2 !!)
# parameters[1, :] = - xi
# parameters[2, :] = sigma
############################################

### these functions take untransformed parameters as arguments

def firstStateGenerator(parameters, size):
    first_state = zeros((size, 2))
    first_state[:, 0] = random.normal(size = size, loc = 520, scale = 10)
    return first_state
def transitionAndWeight(states, y, parameters, t):
    Nx = states.shape[0]
    Ntheta = states.shape[2]
    noise = random.multivariate_normal(size = (Nx, Ntheta), mean = zeros(2), cov = [[1/3, 0.5], [0.5, 1]])
    weights = zeros((Nx, Ntheta))
    code = \
    """
    float innerterm0 = 0.;
    float innerterm1 = 0.;
    float logdensity0 = 0.;
    float logdensity1 = 0.;
    float logsurv0 = 0.;
    float logsigma = 0.;
    float oneoverxi = 0.;
    for (int j = 0; j < Ntheta; j++)
    {
        logsigma = log(parameters(2, j));
        oneoverxi = - 1 / parameters(1, j);
        for(int k = 0; k < Nx; k++)
        {
            states(k, 0, j) = states(k, 0, j) + states(k, 1, j) + parameters(0, j) * noise(k, j, 0);
            states(k, 1, j) =                   states(k, 1, j) + parameters(0, j) * noise(k, j, 1);
            innerterm0 = 1. + parameters(1, j) / parameters(2, j) * (y(0) - states(k, 0, j));
            innerterm1 = 1. + parameters(1, j) / parameters(2, j) * (y(1) - states(k, 0, j));
            logdensity0 = - logsigma - (oneoverxi + 1) * log(innerterm0) - exp(- oneoverxi * log(innerterm0));
            logdensity1 = - logsigma - (oneoverxi + 1) * log(innerterm1) - exp(- oneoverxi * log(innerterm1));
            logsurv0 = - exp(- oneoverxi * log(innerterm0));
            weights(k, j) = logdensity0 + logdensity1 - logsurv0;
        }
    }
    """
    #old surv0, it was wrong: //logsurv0 = log(1 - exp(-exp((-1 / parameters(1, j)) * log(innerterm0))));
    weave.inline(code,['Nx', 'Ntheta', 'states', 'y', 'parameters', 'weights', 'noise'], \
            type_converters=weave.converters.blitz, libraries = ["m"])
    return {"states": states , "weights": weights}


modelx = SSM("Athletics records", xdimension = 2, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setTransitionAndWeight(transitionAndWeight)
modelx.excludedobservations = [17]

###### old functions passed to the old filtering system
#def firststate(xparticles, thetaparticles, t):
#    return xparticles[:, 0, :]
#def CIlow(xparticles, thetaparticles, t):
#    Nx = xparticles.shape[0]
#    Ntheta = xparticles.shape[2]
#    result = zeros((Nx, Ntheta))
#    tempvalue = log(1 - 0.025)
#    for j in range(Ntheta):
#        result[:,j] = xparticles[:, 0, j] - thetaparticles[2, j] / thetaparticles[1, j] * (1 - (- tempvalue)**(+thetaparticles[1, j]))
#    return result
#def CIhigh(xparticles, thetaparticles, t):
#    Nx = xparticles.shape[0]
#    Ntheta = xparticles.shape[2]
#    result = zeros((Nx, Ntheta))
#    tempvalue = log(1 - 0.975)
#    for j in range(Ntheta):
#        result[:, j] = xparticles[:, 0, j] - thetaparticles[2, j] / thetaparticles[1, j] * (1 - (- tempvalue)**(+thetaparticles[1, j]))
#    return result
#def beat1985(xparticles, thetaparticles, t):
#    Nx = xparticles.shape[0]
#    Ntheta = xparticles.shape[2]
#    result = zeros((Nx, Ntheta))
#    for j in range(Ntheta):
#        xioversigma = - thetaparticles[1, j] / thetaparticles[2, j]
#        inner1985 = maximum(0, 1 - xioversigma * (502.62 - xparticles[:, 0, j]))
#        beat1985 = 1 - exp(-(inner1985) ** (+1 / thetaparticles[1, j]))
#        result[:, j] = beat1985 
#    return result
#def beat1993(xparticles, thetaparticles, t):
#    Nx = xparticles.shape[0]
#    Ntheta = xparticles.shape[2]
#    result = zeros((Nx, Ntheta))
#    for j in range(Ntheta):
#        xioversigma = - thetaparticles[1, j] / thetaparticles[2, j]
#        inner1993 = maximum(0, 1 - xioversigma * (486.11 - xparticles[:, 0, j]))
#        beat1993 = 1 - exp(-(inner1993) ** (+1 / thetaparticles[1, j]))
#        result[:, j] = beat1993 
#    return result
#
#modelx.functionals = {"firststate": firststate, "CIlow": CIlow, "CIhigh": CIhigh, \
#        "beat1985": beat1985, "beat1993": beat1993}


