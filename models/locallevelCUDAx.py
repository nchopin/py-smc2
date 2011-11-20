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

transitionKernelTemplate = """
        __global__ void transition(float* y, float* states, float *parameters, float* noise, int* Nx, int* Ntheta, int* statedim)
        {
            int i = threadIdx.x + blockDim.x * blockIdx.x;
            int j = threadIdx.y + blockDim.y * blockIdx.y;
            if (i < *Nx && j < *Ntheta){ 
                int idim = 0;
                int index = idim * (*Ntheta) + i * (*Ntheta * *statedim) + j;
                states[index] = states[index] + sqrt(parameters[0 * (*Ntheta) + j]) * noise[i * (*Ntheta) + j];
                noise[i * (*Ntheta) + j] = -0.9189385 - 0.5 * log(parameters[1 * (*Ntheta) + j]) - 
                0.5 / parameters[1 * (*Ntheta) + j] * (*y - states[index]) * (*y - states[index]);
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
            drv.In(array(Nx, dtype = int32)), \
            drv.In(array(Ntheta, dtype = int32)), \
            drv.In(array(statedim, dtype = int32)), \
            block = (THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, 1), grid = (num_blocks_x, num_blocks_y))
    return {"states": states, "weights": noise}

modelx = SSM(name = "Linear Gaussian model x (using CUDA)", xdimension = 1, ydimension = 1)
modelx.setFirstStateGenerator(firstStateGenerator)
modelx.setObservationGenerator(observationGenerator)
modelx.setTransitionAndWeight(transitionCUDA)
# Values used to generate the synthetic dataset when needed:
# (untransformed parameters)
modelx.parameters = array([0.8, 0.4])
