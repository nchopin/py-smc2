from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, newaxis
from numpy import max as numpymax
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
import pycuda.autoinit
import pycuda.driver as drv
import numpy

from pycuda.compiler import SourceModule

#! /usr/bin/env python
# -*- coding: utf-8 -*-


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
            //term = states(k, 0, j);
            // * states(k, 0, j) / 20)
            term = (double) y(0) - (states(k, 0, j) * states(k, 0, j) / 20);
            weights(k, j) = tempmeasure1 + tempmeasure2 * term * term;
            //weights(k, j) = term;
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

Nx = 4
Ntheta = 3
statedim = 2
random.seed(820)
states = random.normal(size = (Nx, statedim, Ntheta))
y = 1
parameters = ones((2, Ntheta))
random.seed(820)
transitionKernelTemplate = """
        __global__ void transition(float* y, float* states, float *parameters, float* noise, float* t, int* Nx, int* Ntheta, int* statedim)
        {
            int i = threadIdx.x + blockDim.x * blockIdx.x;
            int j = threadIdx.y + blockDim.y * blockIdx.y;
            if (i < *Nx && j < *Ntheta){ 
                float term = 0;
                //int idim = 0;
                int index = 0 * (*Ntheta) + i * (*Ntheta * *statedim) + j;
                //term = y[0] - (states[index] * states[index] / 20.);
                states[index] = 0.5 * states[index] + 25 * (states[index] / (1 + states[index] * states[index])) 
                + 8 * cosf(1.2 * (t[0] - 1)) + sqrt(parameters[0 * (*Ntheta) + j]) * noise[i * (*Ntheta) + j];
                //term = states[index];
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


print "CUDA:"
import cProfile
cProfile.run("""
random.seed(820)
CPUres = transitionAndWeight(states.copy(), y, parameters, 4)
random.seed(820)
cudares = transitionCUDA(states, y, parameters, 4)
""", "prof")
import pstats
p = pstats.Stats('prof')
p.sort_stats('cumulative').print_stats(10)
p.sort_stats('time').print_stats(10)
print "diff"
print numpymax(CPUres["states"] - cudares["states"])
print numpymax(CPUres["weights"] - cudares["weights"])

print "CPU"
print CPUres["states"]
print "cuda"
print cudares["states"]
print "CPU"
print CPUres["weights"]
print "cuda"
print cudares["weights"]


