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
from src.statespacemodel import SSM
import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule

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
#def transitionAndWeight(states, y, parameters, t):
#    code = \
#    """
#    float tempmeasure1;
#    float tempmeasure2;
#    float temptransition;
#    float term;
#    for (int j = 0; j < Ntheta; j++)
#    {
#        tempmeasure1 = -0.9189385 - 0.5 * log(parameters(1, j));
#        tempmeasure2 = -0.5 / parameters(1, j);
#        temptransition = sqrt(parameters(0, j));
#        for(int k = 0; k < Nx; k++)
#        {
#            states(k, 0, j) = 0.5 * states(k, 0, j) + 25 * states(k, 0, j) /
#            (1 + states(k, 0, j) * states(k, 0, j)) + 8 * cos(1.2 * (t(0) - 1))
#            + temptransition * noise(k, j);
#            term = (double) y(0) - (states(k, 0, j) * states(k, 0, j) / 20);
#            weights(k, j) = tempmeasure1 + tempmeasure2 * term * term;
#        }
#    }
#    """
#    y = array([y])
#    t = array([t])
#    Nx = states.shape[0]
#    Ntheta = states.shape[2]
#    weights = zeros((Nx, Ntheta))
#    noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1)
#    weave.inline(code,['Nx', 'Ntheta', 'states', 'y', 't', 'parameters', 'noise', 'weights'], type_converters=weave.converters.blitz, libraries = ["m"])
#    return {"states": states , "weights": weights}
transitionKernelTemplate = """
        __global__ void transition(float* y, float* states, float *parameters, float* noise, float* t, int* Nx, int* Ntheta, int* statedim)
        {
            int i = threadIdx.x + blockDim.x * blockIdx.x;
            int j = threadIdx.y + blockDim.y * blockIdx.y;
            if (i < *Nx && j < *Ntheta){ 
                float term = 0.;
                int index = 0 * (*Ntheta) + i * (*Ntheta * *statedim) + j;
                states[index] = 0.5 * states[index] + 25 * (states[index] / (1 + states[index] * states[index])) 
                + 8 * cosf(1.2 * (t[0] - 1)) + sqrt(parameters[0 * (*Ntheta) + j]) * noise[i * (*Ntheta) + j];
                term = y[0] - (states[index] * states[index] / 20);
                noise[i * (*Ntheta) + j] = -0.9189385 - 0.5 * log(parameters[1 * (*Ntheta) + j]) - 
                0.5 / parameters[1 * (*Ntheta) + j] * term * term;
            }
        }
"""
#comboKernel = comboKernelTemplate % {"Ntheta": Ntheta}
transitionKernel = transitionKernelTemplate
transitionGF = SourceModule(transitionKernel).get_function("transition")
def transitionCUDA(states, y, parameters, t):
    THREADS_PER_BLOCK_X = 16 
    THREADS_PER_BLOCK_Y = 16 
    y = array(y, dtype = float32)
    states = states.astype(float32)
    parameters = parameters.astype(float32)
    Nx, statedim, Ntheta = states.shape
    noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1).astype(float32)
    num_blocks_x = int(math.ceil(Nx / THREADS_PER_BLOCK_X))
    num_blocks_y = int(math.ceil(Ntheta / THREADS_PER_BLOCK_Y))
    transitionGF(drv.In(y), drv.InOut(states), \
            drv.In(parameters), \
            drv.InOut(noise), \
            drv.In(array(t, dtype = float32)), \
            drv.In(array(Nx, dtype = int32)), \
            drv.In(array(Ntheta, dtype = int32)), \
            drv.In(array(statedim, dtype = int32)), \
            block = (THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, 1), grid = (num_blocks_x, num_blocks_y))
    return {"states": states, "weights": noise}


modelx = SSM(name = "Periodic Gaussian model x (using CUDA)", xdimension = 1, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionCUDA)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.parameters = array([10, 1])




